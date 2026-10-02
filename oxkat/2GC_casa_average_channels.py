# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

# Channel-average a field's MS at the end of 2GC. The input MS holds the self-
# calibrated CORRECTED_DATA and the final model predicted into MODEL_DATA; the
# output holds DATA, CORRECTED_DATA and MODEL_DATA averaged down to the requested
# number of channels, with WEIGHT_SPECTRUM, as in PRE_casa_average_to_1k_add_wtspec.py.

import os
import shutil
import sys


def main():
    # Script arguments follow the script path after `-c`
    args = sys.argv[sys.argv.index('-c') + 2:] if '-c' in sys.argv else sys.argv[1:]
    if len(args) < 3:
        sys.exit('Usage: casa -c oxkat/2GC_casa_average_channels.py <input.ms> <output.ms> <nchan_out>')

    input_ms = args[0].rstrip('/')
    output_ms = args[1].rstrip('/')
    nchan_out = int(args[2])

    if not os.path.isdir(input_ms):
        sys.exit(f'Input MS does not exist: {input_ms}')

    tb.open(input_ms+'/SPECTRAL_WINDOW')
    nchan = int(tb.getcol('NUM_CHAN')[0])
    tb.done()

    tb.open(input_ms)
    in_columns = tb.colnames()
    tb.done()

    if nchan < nchan_out:
        sys.exit(f'{input_ms} has {nchan} channels, fewer than the {nchan_out} requested')

    chanbin = nchan // nchan_out
    print(f'\n{input_ms}: {nchan} channels -> {nchan_out} channels (chanbin = {chanbin})\n')
    if nchan % nchan_out != 0:
        print(f'WARNING: {nchan} channels is not a multiple of {nchan_out}; the output may not have exactly {nchan_out}\n')
    if 'MODEL_DATA' not in in_columns:
        print('WARNING: the input MS has no MODEL_DATA column, so the output will have no predicted model\n')

    # An earlier product of this field is replaced
    for old in (output_ms, output_ms+'.flagversions'):
        if os.path.isdir(old):
            print(f'Removing the existing {old}')
            shutil.rmtree(old)

    out_dir = os.path.dirname(output_ms)
    if out_dir != '' and not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    mstransform(vis = input_ms,
                outputvis = output_ms,
                datacolumn = 'all',
                chanaverage = chanbin > 1,
                chanbin = chanbin,
                realmodelcol = True,
                usewtspectrum = True)

    flagmanager(vis = output_ms, mode = 'save', versionname = 'post-average')

    # Check what came out
    tb.open(output_ms+'/SPECTRAL_WINDOW')
    nchan_made = int(tb.getcol('NUM_CHAN')[0])
    tb.done()

    tb.open(output_ms)
    out_columns = tb.colnames()
    tb.done()

    print(f'\n{output_ms}: {nchan_made} channels; columns {[c for c in out_columns if c in ("DATA", "CORRECTED_DATA", "MODEL_DATA", "WEIGHT_SPECTRUM")]}\n')
    if nchan_made != nchan_out:
        print(f'WARNING: {nchan_made} channels were made, not the {nchan_out} requested\n')
    for column in ('CORRECTED_DATA', 'MODEL_DATA'):
        if column in in_columns and column not in out_columns:
            print(f'WARNING: {column} is in the input MS but not in the output\n')

    clearstat()
    clearstat()


if __name__ == '__main__':
    main()
