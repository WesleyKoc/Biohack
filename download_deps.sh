#!/bin/bash
set -e

echo "Downloading glpk.js..."
curl -L "https://cdn.jsdelivr.net/npm/glpk.js/dist/glpk.min.js" -o glpk.min.js
echo "glpk.min.js: $(wc -c < glpk.min.js) bytes"

echo ""
echo "Done. Still needed manually: js/escher-fba.min.js"
echo "See README.md for instructions."
