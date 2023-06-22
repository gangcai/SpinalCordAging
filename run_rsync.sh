#!/bin/bash
from="/home/gangcai/projects/aging/SpinalCord/analysis_SeuratV5_version2/"
to="SpinalCordScripts"
rsync -zarv --include="*/" --include="*.R" --exclude="*" "$from" "$to"
rsync -zarv --include="*/" --include="*.py" --exclude="*" "$from" "$to"
