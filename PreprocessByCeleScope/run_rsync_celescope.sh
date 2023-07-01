#!/bin/bash
from="/home/db/private/XieLab/aging/mouse_spinal_cord/run_celescope1.11.0/"
to="./"
rsync -zarv --include="*/" --include="*.sh" --exclude="*" "$from" "$to"
