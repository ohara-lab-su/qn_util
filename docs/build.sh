#!/usr/bin/env bash
set -euo pipefail
DOCS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$DOCS_DIR/.." && pwd)"
export PYTHONPATH="$ROOT_DIR/src:${PYTHONPATH:-}"
cd "$DOCS_DIR"
rm -rf _build
for lang in en ja; do
  mkdir -p "$lang/api"
  sphinx-apidoc -f --module-first -o "$lang/api" "$ROOT_DIR/src/qn_util"
  sphinx-build -b html -q -D language="$lang" -D root_doc=index -D exclude_patterns="_build,ja/**,en/**" "$lang" "_build/html/$lang" -c "$DOCS_DIR"
done
# Root page is the English landing page, identical to /en/index.html.
cp "_build/html/en/index.html" "_build/html/index.html"
# Copy the English landing page assets references resolve under /en/, so use a redirect instead.
cat > "_build/html/index.html" <<'EOF'
<!doctype html><html lang="en"><head><meta charset="utf-8"><meta http-equiv="refresh" content="0; url=en/index.html"><link rel="canonical" href="en/index.html"><title>Qn-Util</title></head><body><a href="en/index.html">English documentation</a></body></html>
EOF
