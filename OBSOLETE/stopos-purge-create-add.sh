#!/bin/bash

STOPOS_POOL=d8c24f78f9772cbdff54cf62
stopos purge -p $STOPOS_POOL
stopos create -p $STOPOS_POOL
stopos add -p $STOPOS_POOL tokens/*.json

