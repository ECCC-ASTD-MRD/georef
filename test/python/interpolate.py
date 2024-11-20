#!/usr/bin/env python3

import os
import georef

import argparse
import rmn

""" Testing program for libgeoref

This emulates the program test/C/interpolate.c using Python.  This gives us
a very minimal set of C functions and methods to make a ctypes binding for.

"""

def get_args():
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("-i", "--input",  help="Input file", required=True)
    p.add_argument("-o", "--output", help="Output file", required=True)
    p.add_argument("-t", "--truth",  help="Truth data file to compare with")
    p.add_argument("-g", "--grid",   help="Grid file", required=True)
    p.add_argument("-m", "--type",   choices=['N', 'L', 'V'], default='N', help="Interpolation type")
    p.add_argument("-e", "--etiket", help="ETIKET for destination field")
    p.add_argument("-n", "--nomvar", nargs='*', help="List of variable to process")

    return p.parse_args()

args = get_args()

def main():
    interp = {
        'N': georef.RefInterpR.IR_NEAREST,
        'L': georef.RefInterpR.IR_LINEAR,
        'V': georef.RefInterpR.IR_CUBIC
    }[args.type]
    interpolate(args.input, args.output, args.truth, args.grid, args.nomvar, args.etiket, interp)

# TODO: Librmn already has this function which I could expose to Python but
# I wanted to have it immeidately
def copy_metadata(rec, template):
    rec.ni = template.ni
    rec.nj = template.nj
    rec.nk = template.nk

    rec.deet = template.deet
    rec.npas = template.npas

    rec.etiket = template.etiket
    rec.nomvar = template.nomvar
    rec.grtyp = template.grtyp
    rec.typvar = template.typvar

    rec.ip1 = template.ip1
    rec.ip2 = template.ip2
    rec.ip3 = template.ip3

    rec.ig1 = template.ig1
    rec.ig2 = template.ig2
    rec.ig3 = template.ig3
    rec.ig4 = template.ig4

    rec.dateo = template.dateo
    rec.datev = template.datev

def interpolate(input_file: str, output_file: str, truth: str, grid: str, vars: [str], etiket: str, interp_type: int):

    input_file = rmn.fst24_file(filename=input_file)

    output_file = rmn.fst24_file(filename=output_file, options='R/W')

    grid_file = rmn.fst24_file(filename=grid)

    options_1 = georef.TGeoOptions(Interp=georef.RefInterpR.IR_NEAREST, NoData=0.0)


    grid_rec = next(grid_file.new_query(nomvar="GRID"))
    refout = georef.TGeoref(grid_rec.ni, grid_rec.nj, grid_rec.grtyp, grid_rec.ig1, grid_rec.ig2, grid_rec.ig3, grid_rec.ig4, grid_file)
    refout.valid()

    for rec_in in input_file.new_query(nomvar=args.nomvar[0]):
        refin = georef.TGeoref(rec_in.ni, rec_in.nj, rec_in.grtyp, rec_in.ig1, rec_in.ig2, rec_in.ig3, rec_in.ig4, input_file)
        refin.valid()

        georef.Interp(refout, refin, options_1, grid_rec.data, rec_in.data)
        copy_metadata(grid_rec, rec_in)
        output_file.write(grid_rec, rewrite=True)

main()
