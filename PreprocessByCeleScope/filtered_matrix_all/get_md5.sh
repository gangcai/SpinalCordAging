#!/bin/bash
find ./SC* -type f -exec md5sum {} + > matrix_files.md5
echo "completed"
