#!/bin/bash
# Automatic Git backup script

echo "Backing up to GitHub..."

git add .

if [ -z "$1" ]; then
    git commit -m "Auto-backup: $(date '+%Y-%m-%d %H:%M:%S')"
else
    git commit -m "$1"
fi

git push origin main

echo "Backup complete!"
echo "View at: https://github.com/nsolly03/rtc-ppi-inhibitor-discovery"

