#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
DESCRIPTION:  Code to run a n-fold cross validation on the results of the GRASS
              GIS v.surf.bspline and v.surf.idw function. This code is used in
              a tutorial about carrying out n-fold cross validation in GRASS
              GIS (https://tutorials.ecodiv.earth/toc/cross_validation.html.
NOTE:         Code should work on GRASS GIS 7.2 + and should be run from
              within a GRASS session.

@author: pvbreugel add ecodiv dot earth (2016)
"""

# Modules
# -----------------------------------------------------------------------------
import os
import sys
import numpy as np
import uuid
import tempfile
import string
from grass.pygrass.modules import Module
import grass.script as gs
from subprocess import PIPE


# Functions
# -----------------------------------------------------------------------------
def tmpname(prefix):
    """Generate a tmp name which contains prefix
    Store the name in the global list.
    Use only for raster maps.
    """
    tmpf = prefix + str(uuid.uuid4())
    tmpf = string.replace(tmpf, '-', '_')
    return tmpf


def bspline_param(vectormap, depvar):
    """Get output bspline parameter estimates"""
    stfact = Module("v.surf.bspline", flags="ec", input=vectormap,
                    column=depvar, memory=1024, stdout_=PIPE).outputs.stdout
    stepdist = float(stfact.split(':')[-1].strip())
    stfact = Module("v.surf.bspline", flags="c", input=vectormap,
                    ew_step=stepdist, ns_step=stepdist, column=depvar,
                    memory=1024, stdout_=PIPE).outputs.stdout
    stfact = stfact.replace(" ", "")
    stfact = stfact.split("\n")[1:-1]
    stfact = [z.split("|") for z in stfact]
    stfact = [[float(x) for x in y if x != ''] for y in stfact]
    minlambda = min(stfact, key=lambda x: abs(x[1]))[0]
    return(stepdist, minlambda)


def bspline_validation(vector, column, lambda_i, ew_step, ns_step, keep,
                       npartitions=4, method='bilinear', solver='cholesky',
                       maxit=10000, memory=2048):
    """Compute validation statistics (rsme) for bspline extrapolation"""

    # Temporary map
    tmpmap = tmpname("cvbase_")
    Module("g.copy", vector=[vector, tmpmap])

    # Compute rsme over model with all callibration points
    tmpout = tmpname("cvtmp4_")
    Module("v.surf.bspline", input=tmpmap, column=column, ew_step=ew_step,
           ns_step=ns_step, method=method, lambda_i=lambda_i,
           solver=solver, maxit=maxit, memory=memory, raster_output=tmpout)
    Module("v.what.rast", map=tmpmap, raster=tmpout, column="bspline")
    stats = Module("db.select", flags="c", sql="SELECT {},bspline FROM {}".
                   format(column, tmpmap), stdout_=PIPE).outputs.stdout
    stats = stats.replace("\n", "|")[:-1].split("|")
    stats = (np.asarray([float(x) for x in stats], dtype="float").
             reshape(len(stats)/2, 2))
    rsme_all = np.sqrt(np.mean(np.diff(stats, axis=1)**2))
    if keep:
        Module("g.rename", raster=[tmpout, keep])
    else:
        Module("g.remove", type="raster", name=tmpout, flags="f")

    # Run n-fold crossvalidation
    if npartitions > 0:
        Module("v.kcv", map=tmpmap, npartitions=npartitions)
        rsme = []
        for i in range(1, npartitions+1):
            tmp_cal = tmpname("cvtmp_calibrate_")
            tmp_val = tmpname("cvtmp_validate_")
            tmpout1 = tmpname("cvtmp_output_1_")
            tmpout2 = tmpname("cvtmp_output_2_")
            tmpout3 = tmpname("cvtmp_output_3_")
            Module("v.extract", flags="r", input=tmpmap, output=tmp_cal,
                   where="part={}".format(i))
            Module("v.extract", input=tmpmap, where="part={}".format(i),
                   output=tmp_val)
            Module("v.surf.bspline", input=tmp_cal, column=column,
                   ew_step=ew_step, ns_step=ns_step, method=method,
                   lambda_i=lambda_i, solver=solver, maxit=maxit,
                   memory=memory, output=tmpout1, sparse_input=tmp_val)
            Module("v.category", input=tmpout1, output=tmpout2,
                   option="del", cat=-1)
            Module("v.category", input=tmpout2, output=tmpout3, option="add")
            Module("v.db.addtable", map=tmpout3)
            Module("v.db.addcolumn", map=tmpout3,
                   columns=("x double precision, y double precision, "
                            "z double precision"))
            Module("v.to.db", map=tmpout3, option="coor", columns="x,y,z")
            # TODO: need to find out how to use the from_ with Module
            gs.run_command("v.distance", from_=tmpout3, to=tmp_val,
                           upload="to_attr", column="x", to_column=column)
            stats = Module("db.select", flags="c", sql="SELECT x, z FROM {}".
                           format(tmpout3), stdout_=PIPE).outputs.stdout
            stats = stats.replace("\n", "|")[:-1].split("|")
            stats = (np.asarray([float(x) for x in stats], dtype="float").
                     reshape(len(stats)/2, 2))
            rsme.append(np.sqrt(np.mean(np.diff(stats, axis=1)**2)))
            Module("g.remove", type="vector", pattern="cvtmp_*", flags="f")
        Module("g.remove", type="vector", pattern="cvbase_*", flags="f")
        return {'rsme_full': rsme_all, 'rsme_cv_mean': np.asarray(rsme).mean(),
                'rsme_cv_std': np.asarray(rsme).std(), 'rsme_cv': rsme}
    else:
        return {'rsme_full': rsme_all}


def idw_validation(vector, column, keep, npoints=12, power=2, npartitions=10,
                   memory=2048):
    """Compute validation statistics (rsme) for idw extrapolation"""

    # Temporary map
    tmpmap = tmpname("cvbase_")
    Module("g.copy", vector=[vector, tmpmap])

    # Compute rsme over model with all callibration points
    tmpout = tmpname("cvtmp4_")
    Module("v.surf.idw", input=tmpmap, column=column, npoints=npoints,
           power=power, output=tmpout)
    Module("v.what.rast", map=tmpmap, raster=tmpout, column="idw")
    stats = Module("db.select", flags="c", sql="SELECT {},idw FROM {}".
                   format(column, tmpmap), stdout_=PIPE).outputs.stdout
    stats = stats.replace("\n", "|")[:-1].split("|")
    stats = (np.asarray([float(x) for x in stats], dtype="float").
             reshape(len(stats)/2, 2))
    rsme_all = np.sqrt(np.mean(np.diff(stats, axis=1)**2))
    if keep:
        Module("g.rename", raster=[tmpout, keep])
    else:
        Module("g.remove", type="raster", name=tmpout, flags="f")

    # Run n-fold crossvalidation
    if npartitions > 0:
        Module("v.kcv", map=tmpmap, npartitions=npartitions)
        rsme = []
        for i in range(1, npartitions+1):
            tmppnt = tmpname("cvtmp2_")
            tmpspa = tmpname("cvtmp3_")
            tmpout = tmpname("cvtmp4_")
            Module("v.extract", flags="r", input=tmpmap, output=tmppnt,
                   where="part={}".format(i))
            Module("v.extract", input=tmpmap, where="part={}".format(i),
                   output=tmpspa)
            Module("v.surf.idw", input=tmppnt, column=column, npoints=npoints,
                   power=power, output=tmpout)
            Module("v.what.rast", map=tmpspa, raster=tmpout, column="idw")
            stats = Module("db.select", flags="c",
                           sql="SELECT {},idw FROM {}".
                           format(column, tmpspa), stdout_=PIPE).outputs.stdout
            stats = stats.replace("\n", "|")[:-1].split("|")
            stats = (np.asarray([float(x) for x in stats], dtype="float").
                     reshape(len(stats)/2, 2))
            rsme.append(np.sqrt(np.mean(np.diff(stats, axis=1)**2)))
            Module("g.remove", type="all", pattern="cvtmp*", flags="f")
        Module("g.remove", type="vector", pattern="cvbase_*", flags="f")

    # Return output
    return {'rsme_full': rsme_all, 'rsme_cv_mean': np.asarray(rsme).mean(),
            'rsme_cv_std': np.asarray(rsme).std(), 'rsme_cv': rsme}

# Example
# -----------------------------------------------------------------------------

# Determine parameters
stepdist, minlambda = bspline_param(vectormap="households2", depvar="lv")

# Compute evaluation statistics
bspline_stats = bspline_validation(vector="households2", column="lv",
                                   keep="lv2_bspline", lambda_i=minlambda,
                                   ew_step=stepdist, ns_step=stepdist,
                                   npartitions=10, method='bilinear',
                                   solver='cholesky', maxit=10000, memory=2048)
idw_validation(vector="households2", column="lv", npoints=12, power=2,
               npartitions=10, keep="lv2_idw")


































