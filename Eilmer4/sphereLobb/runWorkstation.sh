#!/bin/bash

ls -d */ | parallel --bar --jobs 4 'cd {} && ./run-calculation-in-stages.lua'