#!/bin/bash
################################################################################
# File:        lago-anna.sh
# Description: Defines the environmental variables and modify the $USER .bashrc
# Author:      Hernán Asorey
# Email:       asoreyh@gmail.com
# Date:        2012-2023
# 
# Copyright:   [2012] [The LAGO Collaboration]
# License:     BSD-3-Clause
# See the LICENSE file in the project root for full license information.
################################################################################
VERSION="1.5"
export LAGO_ANNA=${PWD}
export LAGO_ANNA_VERSION=${VERSION}
date=$(date -u)
echo "#
## Changes added by the LAGO ANNA suite on $date
#
export LAGO_ANNA=\"${LAGO_ANNA}\"
export LAGO_ANNA_VERSION=\"${VERSION}\"
export PATH=\"\${LAGO_ANNA}:\$PATH\"
" >> ${HOME}/.bashrc
source ${HOME}/.bashrc