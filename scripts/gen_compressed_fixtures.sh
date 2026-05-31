#!/usr/bin/env bash
#
# Generate tiled-image-compressed FITS fixtures for fits4 integration tests.
#
# fits4 reads the FITS tiled-image compression convention (RICE_1, GZIP_1,
# GZIP_2, PLIO_1, HCOMPRESS_1) in pure Rust. To test that decoder against the
# canonical encoder, we use `fpack` (shipped with cfitsio) to compress known
# source images, then assert in our integration tests that decompressing them
# reproduces the original pixels. `fpack` is a *dev/test* dependency only — the
# fits4 crate itself stays pure Rust.
#
# Source images are the uncompressed NASA samples already in samp/. Outputs are
# written to samp/ as <stem>.<algo>.fits.fz and round-trip-verified with funpack.
#
# Requirements (macOS): brew install cfitsio   (provides fpack + funpack)
#            (conda):    conda install -c conda-forge cfitsio
# Upload requires the gcloud SDK (`gcloud storage` or `gsutil`) with write
# access to the fits4_samples bucket.
#
# Usage:
#   scripts/gen_compressed_fixtures.sh            # generate + verify locally
#   scripts/gen_compressed_fixtures.sh --upload   # also upload fixtures to GCS
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SAMP_DIR="$REPO_ROOT/samp"
BUCKET="gs://fits4_samples"

UPLOAD=0
[[ "${1:-}" == "--upload" ]] && UPLOAD=1

# --- tool checks ------------------------------------------------------------
for tool in fpack funpack; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "error: '$tool' not found on PATH. Install cfitsio:" >&2
    echo "  macOS:  brew install cfitsio" >&2
    echo "  conda:  conda install -c conda-forge cfitsio" >&2
    exit 1
  fi
done

if [[ ! -d "$SAMP_DIR" ]]; then
  echo "error: $SAMP_DIR not found. Fetch the base samples first." >&2
  exit 1
fi

GENERATED=()

# compress <source-stem.fits> <algo-label> <fpack-flags...>
#   produces samp/<stem>.<algo-label>.fits.fz and verifies it round-trips.
compress() {
  local src="$1"; shift
  local label="$1"; shift
  local src_path="$SAMP_DIR/$src"
  local stem="${src%.fits}"
  local out="$SAMP_DIR/${stem}.${label}.fits.fz"

  if [[ ! -f "$src_path" ]]; then
    echo "  skip: $src not present in samp/ — skipping $label" >&2
    return 0
  fi

  rm -f "$out"
  # -O names the output explicitly; remaining args are the algorithm flags.
  fpack "$@" -O "$out" "$src_path"

  # Verify: funpack -> /dev/stdout would re-add a primary; instead unpack to a
  # temp file and confirm it succeeds (byte-exactness vs. source is asserted in
  # the Rust tests, which compare fits4's decode against the original samples).
  local tmp; tmp="$(mktemp -t funpack.XXXXXX).fits"
  if funpack -O "$tmp" "$out" >/dev/null 2>&1; then
    local before after
    before=$(wc -c < "$src_path")
    after=$(wc -c < "$out")
    printf "  ok  %-40s %8s -> %8s bytes\n" "$(basename "$out")" "$before" "$after"
    GENERATED+=("$(basename "$out")")
  else
    echo "  FAIL: funpack could not decompress $out" >&2
    rm -f "$tmp"
    return 1
  fi
  rm -f "$tmp"
}

echo "Generating compressed fixtures in $SAMP_DIR ..."

# --- integer images (RICE / GZIP / PLIO) ------------------------------------
# EUVEngc4151imgx.fits has I16 IMAGE extensions.
compress EUVEngc4151imgx.fits     rice   -r
compress EUVEngc4151imgx.fits     gzip1  -g
compress EUVEngc4151imgx.fits     gzip2  -g2
# PLIO_1 targets integer (mask-like) data, <= 24-bit; fine for I16.
compress EUVEngc4151imgx.fits     plio   -p
# Force square tiles to exercise 2-D reassembly with edge-truncated tiles.
compress EUVEngc4151imgx.fits     rice_t100   -r -t 100,100

# FGSf64y0106m_a1f.fits is an I32 image.
compress FGSf64y0106m_a1f.fits    rice   -r
compress FGSf64y0106m_a1f.fits    gzip1  -g

# --- float images (quantize; with and without subtractive dither) -----------
# Note: RICE/HCOMPRESS need quantization for floats. `-q <lvl>` quantizes with
# subtractive dithering; `-q0 <lvl>` (no space) quantizes with NO dithering
# (ZQUANTIZ=NO_DITHER); lossless float (`-q 0`, level zero) is GZIP-only.
# FOCx38i0101t_c0f.fits is F32 1024x1024.
compress FOCx38i0101t_c0f.fits    rice_dith     -r -q 16     # quantize + SUBTRACTIVE_DITHER_1
compress FOCx38i0101t_c0f.fits    rice_dith2    -r -qz5 16   # SUBTRACTIVE_DITHER_2, fixed seed 5
compress FOCx38i0101t_c0f.fits    rice_nodith   -r -q0 16    # quantize, NO_DITHER
compress FOCx38i0101t_c0f.fits    gzip_lossless -g -q 0      # lossless float (no quantization)
compress FOCx38i0101t_c0f.fits    hcomp         -h

# WFPC2u5780205r_c0fx.fits is an F32 200x200x4 cube (multi-axis tiling).
compress WFPC2u5780205r_c0fx.fits rice_dith     -r -q 16
compress WFPC2u5780205r_c0fx.fits rice_nodith   -r -q0 16

echo
echo "Generated ${#GENERATED[@]} fixtures:"
printf '  %s\n' "${GENERATED[@]}"

# --- optional upload --------------------------------------------------------
if [[ "$UPLOAD" -eq 1 ]]; then
  if command -v gcloud >/dev/null 2>&1; then
    GSCMD=(gcloud storage cp)
  elif command -v gsutil >/dev/null 2>&1; then
    GSCMD=(gsutil cp)
  else
    echo "error: --upload requested but neither 'gcloud' nor 'gsutil' is installed." >&2
    exit 1
  fi
  echo
  echo "Uploading to $BUCKET ..."
  for f in "${GENERATED[@]}"; do
    "${GSCMD[@]}" "$SAMP_DIR/$f" "$BUCKET/$f"
    echo "  uploaded $f"
  done
fi

echo
echo "Next steps:"
echo "  1. Add these filenames to the download loop in .github/workflows/ci.yml"
echo "     (and bump the cache key, e.g. fits4-samples-v1 -> v2)."
echo "  2. Add assertions in tests/sample_files.rs using hdu.as_compressed_image()."
if [[ "$UPLOAD" -eq 0 ]]; then
  echo "  3. Re-run with --upload to push fixtures to $BUCKET."
fi
