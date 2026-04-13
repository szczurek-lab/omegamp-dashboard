#!/bin/bash
# Build and optionally password-protect the OmegAMP dashboard.
#
# Usage:
#   ./deploy.sh                    # build unprotected
#   ./deploy.sh --password SECRET  # build + encrypt with password
#
# Requirements: python3, modlamp, node/npm (for npx staticrypt)

set -euo pipefail
cd "$(dirname "$0")"

PASSWORD=""
REMEMBER_DAYS=30

while [[ $# -gt 0 ]]; do
    case $1 in
        --password|-p) PASSWORD="$2"; shift 2 ;;
        --remember)    REMEMBER_DAYS="$2"; shift 2 ;;
        *)             echo "Unknown option: $1"; exit 1 ;;
    esac
done

echo "=== Building dashboard ==="
python3 build_dashboard.py --data-dir data --template dashboard_template.html --output docs/index.html

if [[ -n "$PASSWORD" ]]; then
    echo ""
    echo "=== Encrypting with password ==="
    # Keep unencrypted copy for local development
    cp docs/index.html docs/index_unprotected.html

    npx staticrypt docs/index.html \
        -p "$PASSWORD" \
        --remember "$REMEMBER_DAYS" \
        --short \
        --template-title "OmegAMP Dashboard" \
        --template-instructions "Enter the password shared by the corresponding author." \
        --template-button "Open Dashboard" \
        --template-placeholder "Password" \
        --template-color-primary "#22c55e" \
        --template-color-secondary "#292524" \
        --template-error "Wrong password. Contact the corresponding author." \
        --template-remember "Remember me for ${REMEMBER_DAYS} days" \
        -d docs_enc \
        -c false

    mv docs_enc/index.html docs/index.html
    rmdir docs_enc
    echo "✓ Encrypted. Password: $PASSWORD"
    echo "  Unencrypted copy: docs/index_unprotected.html"
else
    # Remove stale encrypted files
    rm -f docs/index_unprotected.html
    echo "✓ No password set (unprotected)"
fi

echo ""
echo "=== Done ==="
ls -lh docs/index.html
echo ""
echo "To deploy: git add -A && git commit -m 'update dashboard' && git push"
