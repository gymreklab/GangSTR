#!/bin/sh

# Create a new git tag before running this script 
# git tag -a vX.Y.Z -m"message"
die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}
./reconf 
./configure || die "./configure failed"

GITID=$(./config/git-version-gen .tarball-version) || die "Failed to get GIT ID"
TARBALL=GangSTR-${GITID}.tar.gz
rm -f "$TARBALL"

make distcheck || die "make distcheck failed"
[ -e "$TARBALL" ] || die "can't find Tarball file '$TARBALL' after 'make distcheck'"

echo
echo "Version $GITID is ready for distribution"
echo ""
echo "Upload the following file to GitHub:"
echo "  $TARBALL"
echo ""
echo ""
