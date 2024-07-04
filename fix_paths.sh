#!/bin/bash

# Makes sure Snakemake only writes to the current directory
# By default, Snakemake uses cache paths determined by the deprecated `appdirs` package, which ends up being in the
# current user's home dir. Problem is, this doesn't work in a cluster environment (home dir is RO).
# It can be overridden to a local directory as shown below:

mkdir -p .tmp/local
export XDG_DATA_HOME=$(pwd)/.tmp/local
mkdir -p .tmp/config
export XDG_CONFIG_HOME=$(pwd)/.tmp/config
mkdir -p .tmp/cache
export XDG_CACHE_HOME=$(pwd)/.tmp/cache
mkdir -p .tmp/state
export XDG_STATE_HOME=$(pwd)/.tmp/state
