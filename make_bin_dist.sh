#!/bin/sh

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

GITID=$(./config/git-version-gen .tarball-version) || die "Failed to get GIT ID"
TARBALL=GangSTR-${GITID}.tar.gz
[ -e "$TARBALL" ] || die "can't find source code tarball file '$TARBALL' (did you run 'make_src_dist.sh'?)"

MACHINE=$(uname -m) || die "uname -m failed"
KERNEL=$(uname -s) || die "uname -s failed"

# new temporary directory, into which we'll extract the files
BUILD_ROOT_DIR=$(mktemp -d build.XXXXXXX) || die "mktemp failed"
BUILD_ROOT_DIR="$PWD/$BUILD_ROOT_DIR"

# the source code directory, extracted from the Tarball
BUILD_SRC_DIR="$BUILD_ROOT_DIR/GangSTR-${GITID}"

# The target prefix directory for "./configure"
GANGSTR_BIN_PREFIX="GangSTR-bin-${KERNEL}-${MACHINE}-${GITID}"
BUILD_BIN_DIR="$BUILD_ROOT_DIR/$GANGSTR_BIN_PREFIX"

# Extract the source code files
tar -xzf "$TARBALL" -C "$BUILD_ROOT_DIR" || die "tar -xf ../$TARBALL failed"

# Build GangSTR, with special static options
cd "$BUILD_SRC_DIR" || die "cd $BUILD_SRC_DIR failed"
./configure --enable-all-static --prefix "$BUILD_BIN_DIR" || die "Configure failed, check logs in $BUILD_SRC_DIR"
make || die "make failed, check logs in $BUILD_SRC_DIR"
make install || die "make install failed, check logs in $BUILD_SRC_DIR"
## The compiled binary files are now in "$BUILD_BIN_DIR"

# Strip the binaries to reduce their size
find "$BUILD_BIN_DIR/bin" -type f | xargs -L 1 strip
# No need to distribute the M4 macroes
rm -rf "$BUILD_BIN_DIR/share/aclocal"

# Pack the binaries for a distribution
cd "$BUILD_ROOT_DIR" || die "cd $BUILD_ROOT_DIR failed"
BINTARBALL="../${GANGSTR_BIN_PREFIX}.tar.gz"
rm -f "$BINTARBALL"
tar -czf "$BINTARBALL" "$GANGSTR_BIN_PREFIX" || die "tar -czf failed"

## Back to the source code directory
cd .. || die "cd .. failed"
rm -r "$BUILD_ROOT_DIR" || die "cleanup failed (rm -r $BUILD_ROOT_DIR)"

echo
echo "Binary Distribution for GangSTR version $GITID for $KERNEL-$MACHINE is ready"
echo ""
echo "Upload the following file to github/releases:"
echo "  $GANGSTR_BIN_PREFIX.tar.gz"
echo ""
