{
  description = "sgffp — SnapGene File Format Parser";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            python312
            uv
            ruff
          ];

          shellHook = ''
            echo "sgffp dev shell — python $(python3 --version 2>&1 | cut -d' ' -f2), uv $(uv --version 2>&1 | cut -d' ' -f2)"
          '';
        };
      }
    );
}
