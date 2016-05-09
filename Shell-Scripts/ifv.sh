#!/bin/bash
# Script to copy over all necessary files for simulation of system in a vacuum
# Edit files copied in as necessary

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for file in wiggle.mdp NaPore.top em.mdp; do
	cp $DIR/Files/$file $PWD

done

echo 'hello'

