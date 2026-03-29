# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import os
import sys


def main():
    # CASA passes script args after `-c`, so require input and output MS paths.
    if len(sys.argv) < 3:
        sys.exit('Usage: casa -c tools/casa_split_corrected_ms.py <input.ms> <output.ms> [field]')

    input_ms = sys.argv[1]
    output_ms = sys.argv[2]
    field = sys.argv[3] if len(sys.argv) > 3 else ''

    if not os.path.isdir(input_ms):
        sys.exit(f'Input MS does not exist: {input_ms}')

    if os.path.isdir(output_ms):
        sys.exit(f'Output MS already exists: {output_ms}')

    # Split CORRECTED_DATA into a standalone MS for downstream calibration/imaging.
    kwargs = {
        'vis': input_ms,
        'outputvis': output_ms,
        'usewtspectrum': True,
        'realmodelcol': True,
        'datacolumn': 'corrected',
    }

    if field.strip() != '':
        kwargs['field'] = field

    print(f'Splitting CORRECTED_DATA from {input_ms} -> {output_ms}')
    if field.strip() != '':
        print(f'Field selection: {field}')

    mstransform(**kwargs)

    # Save a checkpoint flag version on the newly created MS.
    flagmanager(vis=output_ms,
                mode='save',
                versionname='post-split-corrected')

    print(f'Done: {output_ms}')


if __name__ == '__main__':
    main()
