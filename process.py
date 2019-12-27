#coding: utf-8

import numpy as np
import pandas as pd
import os
from datetime import datetime, timedelta

data_dir = 'E4.5'
methy_point_file_path = 'methypoint.tsv'
df_methyinfo = pd.read_csv(methy_point_file_path, sep='\t', dtype=str)
chr_list = df_methyinfo['chrom'].unique()
df_methyinfo[['region start', 'region end']] = df_methyinfo[['region start', 'region end']].applymap(np.int)
df_methyinfo = df_methyinfo.sort_values(by=['chrom', 'region start', 'region end'], axis=0)
df_methyinfo.reset_index(drop=True, inplace=True)
df_output = df_methyinfo.loc[:, ['chrom', 'region start', 'region end', 'region size']].copy()
range_num = df_output.shape[0]

print('file_num = ', len(os.listdir(data_dir)))
start_time = datetime.now()
print('start time = ', start_time.strftime('%Y-%m-%d %H:%M:%S'))

cnt = 0
for file_name in os.listdir(data_dir):
	file_path = os.path.join(data_dir, file_name)
	cell_name = '_'.join(file_name.split('_')[-4:])[:-4]
	cnt += 1
	print(cnt, ' ', cell_name)
	df_methy = pd.read_csv(file_path, sep=',', header=None, names=['chrom', 'region start', 'region end', 'fraction', 'a', 'b'], dtype=str)
	df_methy['fraction'] = df_methy['fraction'].apply(np.float)
	df_methy[['region start', 'region end']] = df_methy[['region start', 'region end']].applymap(np.int)
	df_methy = df_methy.sort_values(by=['chrom', 'region start', 'region end'], axis=0)
	df_methy.reset_index(drop=True, inplace=True)
	#assert (df_methy['region start'] == df_methy['region end']).sum() == df_methy.shape[0]
	cpg_num_label = cell_name + '_cpgnum'
	cpg_fraction_label = cell_name + '_fraction'
	initial_data = np.zeros((range_num, 2), dtype=np.float)
	df_data = pd.DataFrame(initial_data, columns=[cpg_num_label, cpg_fraction_label])
	df_data[cpg_num_label] = df_data[cpg_num_label].apply(np.int)

	info_idx = 0	
	def get_data(item_series):
		global info_idx
		while info_idx < range_num:
			start = df_output.loc[info_idx, 'region start']
			end = df_output.loc[info_idx, 'region end']
			point = item_series['region start']
			if item_series['chrom'] == df_output.loc[info_idx, 'chrom']:
				if point >= start and point <= end:
					df_data.loc[info_idx, cpg_num_label] += 1
					df_data.loc[info_idx, cpg_fraction_label] += item_series['fraction']
					break
				else:
					if point < start:  #very important
						return
					info_idx += 1
			else:
				if item_series['chrom'] not in chr_list:
					return
				info_idx += 1
		return

	df_methy.apply(get_data, axis=1)
	df_data[cpg_fraction_label] = df_data.apply(lambda x: 0 if x[cpg_num_label] == 0 else x[cpg_fraction_label] / x[cpg_num_label], axis=1)
	df_output = pd.concat([df_output, df_data], axis=1)

end_time = datetime.now()
print('end time = ', end_time.strftime('%Y-%m-%d %H:%M:%S'))
print('time duration = ', int((end_time - start_time).total_seconds() / 60))

output_file_path = data_dir + '_methy_fraction.csv' 
df_output.to_csv(output_file_path, sep=',', index=False)
