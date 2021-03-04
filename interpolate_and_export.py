import pandas as pd
import os
import numpy as np
import gridpp

def interpolate_and_save(
    var_name,
    product,
    date_str,
    ds_cf,
    ds_pf,
    in_grid,
    out_points,
    data_export_dir,
    data_BW,
):
    file_path = os.path.join(data_export_dir, f'{var_name}_{product}_{date_str}_bilinear.csv')

    print('File save location: ' + file_path)

    df_out = pd.DataFrame(
        columns = tuple([
            'locNo', 
            'date',
            'step_fwrd', 
            f'{var_name}_ctrl'
        ] + [
            f'{var_name}_ptrb_{num:02}'
            for num in ds_pf.get_index('number')
        ])
    )

    for step in ds_cf.get_index('step'):
        print('Time step: ' + str(step))
        # NB: Is this the best way to deal with missings from on-land coordinates? 
        # NB: Axes of lat/lon are reversed between gridpp and xarray.
        print('Filling')
        gridpp.fill_missing(np.transpose(ds_cf.variables[var_name][step.days - 1,:,:].data)) 
        print('Done')

        cf_values = gridpp.bilinear(
            in_grid,
            out_points,
            gridpp.fill_missing(np.transpose(ds_cf.variables[var_name][step.days - 1,:,:].data)) 
        )
        pf_values = np.empty((len(data_BW), len(ds_pf.get_index('number'))), dtype=float)

        for num in ds_pf.get_index('number'):
            pf_values[:, num - 1] = gridpp.bilinear(
                in_grid, 
                out_points, 
                gridpp.fill_missing(np.transpose(ds_pf.variables[var_name][num - 1, step.days - 1,:,:].data))
            )
        for i in range(len(data_BW)):
            row_dict = {
                'locNo' : data_BW['localityNo'][i],
                'step_fwrd' : step.days,
                'date' : date_str,
                f'{var_name}_ctrl' : cf_values[i],
            }
            for num in ds_pf.get_index('number'):
                row_dict[f'{var_name}_ptrb_{num:02}'] = pf_values[i, num - 1]

            df_out = df_out.append(pd.Series(row_dict), ignore_index = True)

    df_out.to_csv(file_path)
