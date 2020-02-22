import pandas as pd

df = pd.read_pickle("./pickled_df.pkl")

kernel_filters = ['CM_HIGHER', 'Cumulants', 'CM', 'BGK']
bc_order_filters = ['1st_order_bc', '2nd_order_bc']

size_filters_dict = {'small': 30,
                     'medium': 60,
                     'large': 120
                     }

Pr_filters = [10, 100, 1000]
is_3D = True

if is_3D:
    lattice_size_dict = {'small': '1000x150x3',
                         'medium': '2000x300x3',
                         'large': '4000x600x3'
                         }
    velocity_set = 'D3Q27Q27'
else:
    lattice_size_dict = {'small': '1000x150',
                         'medium': '2000x300',
                         'large': '4000x600'
                         }
    velocity_set = 'D2Q9Q9'

new_df = pd.DataFrame()
legend_df = pd.DataFrame()

for Pr_filter in Pr_filters:
    for size_filter_description in size_filters_dict:
        case_id_filter = 'Pr' + str(Pr_filter) + '_' + size_filter_description

        pre_filtered_df = df[
            (df['D'] == size_filters_dict[size_filter_description]) &
            (df['Pr'] == Pr_filter) &
            (df['is3D'] == is_3D)
            ]

        legend_tmp_df = pd.DataFrame({
            'Case_id': case_id_filter,
            'Lattice Size': lattice_size_dict[size_filter_description],
            'Velocity set': velocity_set,
            'Blockage Ratio': '1/5',
            'D': size_filters_dict[size_filter_description],
            'U': pre_filtered_df['U'].iloc[0],
            'Pr': pre_filtered_df['Pr'].iloc[0],
            'Re': pre_filtered_df['Re'].iloc[0],
            'nu': pre_filtered_df['nu'].iloc[0],
            'k': pre_filtered_df['nu'].iloc[0] / pre_filtered_df['Pr'].iloc[0]
        }, index=[0])

        tmp_df = pd.DataFrame({
            'Case_id': case_id_filter,
            'Nu_' + kernel_filters[0] + '_' + bc_order_filters[0]: None,
            'Nu_' + kernel_filters[0] + '_' + bc_order_filters[1]: None,
            'Nu_' + kernel_filters[1] + '_' + bc_order_filters[0]: None,
            'Nu_' + kernel_filters[1] + '_' + bc_order_filters[1]: None,
            'Nu_FEM': None
        }, index=[0])

        for kernel_filter in kernel_filters:
            for bc_order_filter in bc_order_filters:
                filtered_df = df[
                    (df['Collision_Kernel'] == kernel_filter) &
                    (df['BC_order'] == bc_order_filter) &
                    (df['D'] == size_filters_dict[size_filter_description]) &
                    (df['Pr'] == Pr_filter) &
                    (df['is3D'] == is_3D)
                    ]

                filter_col_name = 'Nu_' + kernel_filter + '_' + bc_order_filter
                tmp_df = tmp_df.assign(**{
                    'Case_id': [case_id_filter],
                    filter_col_name: "%.2f" % ((filtered_df['Nu_avg_source'] + filtered_df['Nu_avg_outlet'])/2.),
                    'Nu_FEM': filtered_df['Nu_FEM']
                })

        new_df = new_df.append(tmp_df)
        legend_df = legend_df.append(legend_tmp_df)

print(new_df)
print(legend_df)

with pd.ExcelWriter(f"LBM_validation_HotKarman_Benchmark{('_2D', '_3D')[is_3D]}.xlsx") as writer:  # doctest: +SKIP
    # dfObj = dfObj.sort_values(by=['is3D', 'Collision_Kernel', 'D', 'BC_order'])
    new_df.to_excel(writer, sheet_name='EnhancedTable', index=False)
    legend_df.to_excel(writer, sheet_name='EnhancedTable', startrow=len(new_df) + 2, index=False)

    # df_Pr10 = df.loc[df['Pr'] == 10]
    df_Pr10 = df.loc[(df['Pr'] == 10) & (df['is3D'] == is_3D)]
    df_Pr10.to_excel(writer, sheet_name='Pr10raw', index=False)
    # df_Pr100 = df.loc[df['Pr'] == 100]
    df_Pr100 = df.loc[(df['Pr'] == 100) & (df['is3D'] == is_3D)]
    df_Pr100.to_excel(writer, sheet_name='Pr100raw', index=False)

    df_Pr1000 = df.loc[(df['Pr'] == 1000) & (df['is3D'] == is_3D)]
    df_Pr1000.to_excel(writer, sheet_name='Pr1000raw', index=False)

print("bye")
