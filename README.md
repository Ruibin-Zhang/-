#Midification
1. Only function_recovery is used.
2. The elevators for different floors is considered for high buildings.
3. The influence of different earthquake intensites on the utility disruption is considered.
4. The resilience is quantified based on R index.

#Pelicun-treads使用流程
	
	1. Pelicun_RID.ipynb #获得维修费用/时间（FEMA P58）
	• 数据准备
		○ 原始数据：EDPs.xlsx
		○ Matlab处理：Monte carlo procedure - SAUSAGE-R1.m
		○ 处理结果：Demand-Sample.xlsx
		○ CMP_marginals.csv
	• 数据处理
		○ Demand模块：添加EDP单位和相应的损失指标
		○ CMP模块：确定性能组和易损性信息
		○ Damage模块：计算各性能组损伤信息
		○ Loss模块：计算维修费用和时间
	• 数据输出
		○ FG_damage.csv 本次选取的易损性数据
		○ FG_cost-time.csv 本次选取的易损性数据对于的修复时间和费用
		○ DMG_sample.csv 各性能组损伤样本
		○ loss_sample.csv 损失样本
		○ DL_summary.csv 损失统计信息
		
	2. Preprocessing.ipynb #数据前处理
	• 数据准备
		○ input_parameters.json 建筑信息-input_parameters.xlsx
	• 数据输出
		○ input_parameters-R.json 修改后的建筑信息
		○ DMG_sample-R.csv 调整表头，用于计算修复费用和时间
		○ DMG_sample-DamageStateWeights.csv 将各损伤状态按照所属极限状态（LS）归类
		○ DV_rec_time.csv 修复时间样本
		○ loss_sample_dmg.csv 损失信息按照易损性组分类
		○ loss_sample_loc.csv 损失信息安照楼层分类
	恢复轨迹程序（数据输出文件夹Output_treads）
	
 	3. Treads_RC
	• 数据准备
		○ Repair_Class_Table.csv 维修等级信息
	• 数据输出
		○ RC_component.csv 将损伤状态与维修等级一一对应
		
	4. Treads_RT
	• 数据准备
		○ IF_delays_input.csv 阻碍因素相关数据
	• 数据输出
		○ RT_RSeq_FR.csv 基于功能恢复的每个修复路径的恢复时间
    
