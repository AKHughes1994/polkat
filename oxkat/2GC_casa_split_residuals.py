# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

# Split a field's residuals into the final MS at the end of 2GC. In the input MS the
# final model has been subtracted from CORRECTED_DATA, so CORRECTED_DATA holds the
# residuals. The output holds them as its DATA column, with WEIGHT_SPECTRUM, and no
# CORRECTED_DATA or MODEL_DATA.

import os
import shutil
import sys


def main():
    # Script arguments follow the script path after `-c`
    args = sys.argv[sys.argv.index('-c') + 2:] if '-c' in sys.argv else sys.argv[1:]
    if len(args) < 2:
        sys.exit('Usage: casa -c oxkat/2GC_casa_split_residuals.py <input.ms> <output.ms>')

    input_ms = args[0].rstrip('/')
    output_ms = args[1].rstrip('/')

    if not os.path.isdir(input_ms):
        sys.exit(f'Input MS does not exist: {input_ms}')

    tb.open(input_ms)
    in_columns = tb.colnames()
    tb.done()

    if 'CORRECTED_DATA' not in in_columns:
        sys.exit(f'{input_ms} has no CORRECTED_DATA column')

    # An earlier product of this field is replaced
    for old in (output_ms, output_ms+'.flagversions'):
        if os.path.isdir(old):
            print(f'Removing the existing {old}')
            shutil.rmtree(old)

    out_dir = os.path.dirname(output_ms)
    if out_dir != '' and not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    print(f'\nSplitting the residuals (CORRECTED_DATA) of {input_ms} -> {output_ms}\n')

    mstransform(vis = input_ms,
                outputvis = output_ms,
                datacolumn = 'corrected',
                usewtspectrum = True)

    flagmanager(vis = output_ms, mode = 'save', versionname = 'post-split-residuals')

    # Check what came out
    tb.open(output_ms)
    out_columns = tb.colnames()
    tb.done()

    print(f'\n{output_ms}: columns {[c for c in out_columns if c in ("DATA", "CORRECTED_DATA", "MODEL_DATA", "WEIGHT_SPECTRUM")]}\n')
    if 'DATA' not in out_columns:
        sys.exit(f'{output_ms} has no DATA column')
    for column in ('CORRECTED_DATA', 'MODEL_DATA'):
        if column in out_columns:
            print(f'WARNING: {column} is in the output, which should hold only the residual DATA\n')

    clearstat()
    clearstat()


if __name__ == '__main__':
    main()
