rm -rf AlphaPhase
mkdir AlphaPhase

# Assumes that the program and manual have both been built.

# To build the program run:
# NOTE: Binaries should be moved to the "binaries" folder and uploaded to bitbucket after builds.
#cmake . ; make

# to build the manual using Sphinx:
( cd alphaphase-doc ; make latexpdf )

cp -r example AlphaPhase

# Copy in the documentation.
cp AlphaPhase-doc/build/latex/AlphaPhase.pdf AlphaPhase/AlphaPhase.pdf

if [ $? != 0 ]; then                   # last command: echo
    echo "The manual needs to be built." # last command: [
    exit 1
fi

# Copy in the binaries
cp binaries/* AlphaPhase

# Create a version file

version=`git describe --tags --abbrev=0`
commit=`git rev-parse --short HEAD`

echo Version: $version > AlphaPhase/version.txt
echo Commit: $commit >> AlphaPhase/version.txt

cp MIT_License.txt AlphaPhase

zip -r AlphaPhase.zip AlphaPhase
