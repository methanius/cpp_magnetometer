{ lib
, stdenv
, cmake
, fetchurl
, openblas
}:

stdenv.mkDerivation rec {
    pname = "armadillo";
    version = "14.0.3";

    src = fetchurl {
        url = "https://sourceforge.net/projects/arma/files/armadillo-14.0.3.tar.xz";
        sha256 = "69YhXusB7kEv7QeMip9/h9Th9hh+vNwbwJ9GCVpPQAM=";

    };

    nativeBuildInputs = [ cmake ];
}
