{
    description = "Nixified master's project.";

    inputs = {
        nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
        flake-utils.url = "github:numtide/flake-utils";
    };

    outputs = {self, nixpkgs, flake-utils,...}@inputs: flake-utils.lib.eachSystem [
    "x86_64-linux"
    ] (system: let
        pkgs = import nixpkgs {
            inherit system;
            };
    in {
        devShells.default = pkgs.mkShell rec {
            name = "Magnetometer master's project";
            packages = with pkgs; [
                gcc
                cmake
                (callPackage ./armadillo.nix {})
            ];
        };

    });
}
