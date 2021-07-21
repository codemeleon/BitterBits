#!/usr/bin/env bash

# Rust and addons
# add `~/.cargo/bin` in .bashrc file
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install du-dust
