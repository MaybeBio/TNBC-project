# 原始数据
data = {
    'BT549': {
        'valid_interaction': 427115856,
        'valid_interaction_rmdup': 310590379,
        'trans_interaction': 35223386,
        'cis_interaction': 275366993,
        'cis_shortRange': 124604747,
        'cis_longRange': 150762246
    },
    'HCC70': {
        'valid_interaction': 587876634,
        'valid_interaction_rmdup': 439133960,
        'trans_interaction': 46475836,
        'cis_interaction': 392658124,
        'cis_shortRange': 129161644,
        'cis_longRange': 263496480
    },
    'HMEC': {
        'valid_interaction': 492430664,
        'valid_interaction_rmdup': 375991517,
        'trans_interaction': 27940991,
        'cis_interaction': 348050526,
        'cis_shortRange': 163813668,
        'cis_longRange': 184236858
    },
    'MB231': {
        'valid_interaction': 235793172,
        'valid_interaction_rmdup': 200168329,
        'trans_interaction': 20676310,
        'cis_interaction': 179492019,
        'cis_shortRange': 50110554,
        'cis_longRange': 129381465
    }
}

# 计算比例并添加新列
for sample in data:
    # 计算比例
    ratios = {
        'valid_interaction_rmdup_ratio': data[sample]['valid_interaction_rmdup'] / data[sample]['valid_interaction'],
        'trans_interaction_ratio': data[sample]['trans_interaction'] / data[sample]['valid_interaction_rmdup'],
        'cis_interaction_ratio': data[sample]['cis_interaction'] / data[sample]['valid_interaction_rmdup'],
        'cis_shortRange_ratio': data[sample]['cis_shortRange'] / data[sample]['cis_interaction'],
        'cis_longRange_ratio': data[sample]['cis_longRange'] / data[sample]['cis_interaction']
    }
    
    # 添加新列
    data[sample].update(ratios)

# 输出结果
for sample, values in data.items():
    print(f"Sample: {sample}")
    for key, value in values.items():
        print(f"{key}: {value:.6f}")
    print()
    
    
    
