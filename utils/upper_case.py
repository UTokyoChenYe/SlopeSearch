patterns = [
    'rrrrrrrrr', 'rrrrrryrr', 'rrrrrryry', 'rrrrrryyr',
    'rrrrrryyy', 'rrryryrrr', 'rrryryyrr', 'rrryryyry',
    'rrryryyyr', 'rrryryyyy', 'rryrrrrrr', 'rryrrryrr',
    'rryrrryry', 'rryrrryyr', 'rryrryrrr', 'rryrryrry',
    'rryrryryr', 'rryrryyrr', 'rryrryyry', 'rryrryyyr',
    'rryrryyyy', 'rryryryyr', 'rryryryyy', 'rryyryrrr',
    'rryyryyrr', 'rryyryyry', 'rryyryyyr', 'rryyryyyy',
    'rryyyryrr', 'rryyyyrrr', 'rryyyyyrr', 'rryyyyyyr',
    'ryrrrryrr', 'ryrrrryry', 'ryrrrryyr', 'ryrrrryyy',
    'ryrrryyyr', 'ryrrryyyy', 'ryrryryyr', 'ryrryryyy',
    'ryryrryrr', 'ryryrryry', 'ryryrryyr', 'ryryrryyy',
    'ryryryrrr', 'ryryryrry', 'ryryryryr', 'ryyrrrrrr',
    'ryyrrrrry', 'ryyrrryrr', 'ryyrrryry', 'ryyrrryyr',
    'ryyrrryyy', 'ryyrryrrr', 'ryyrryrry', 'ryyrryryr',
    'ryyrryyrr', 'ryyrryyry', 'ryyrryyyr', 'ryyrryyyy',
    'ryyryryrr', 'ryyryryry', 'ryyryryyr', 'ryyryryyy',
    'ryyyrrrrr', 'ryyyrrrry', 'ryyyrryrr', 'ryyyrryry',
    'ryyyrryyr', 'ryyyrryyy', 'ryyyryryr', 'ryyyryyrr',
    'ryyyryyry', 'ryyyryyyr', 'ryyyryyyy', 'ryyyyryrr',
    'ryyyyryry', 'ryyyyryyr', 'ryyyyryyy', 'ryyyyyyyr',
    'yrrrrryrr', 'yrrrrryry', 'yrrrrryyr', 'yrrrrryyy',
    'yrryryrrr', 'yryrrrrrr', 'yryrrryrr', 'yryrrryry',
    'yryrrryyr', 'yryryryyr', 'yryryryyy', 'yryyryrrr',
    'yryyryyrr', 'yryyryyry', 'yryyryyyr', 'yryyryyyy',
    'yryyyryrr', 'yryyyyrrr', 'yryyyyyrr', 'yryyyyyyr',
    'yyrrrryrr', 'yyrrrryry', 'yyrrrryyr', 'yyrrrryyy',
    'yyrryryyr', 'yyrryryyy', 'yyryrryrr', 'yyryrryry',
    'yyryrryyr', 'yyryrryyy', 'yyryryrrr', 'yyyrrryrr',
    'yyyrrryry', 'yyyrrryyr', 'yyyrrryyy', 'yyyryryyr',
    'yyyryryyy', 'yyyyrrrrr', 'yyyyrryrr', 'yyyyrryry',
    'yyyyrryyr', 'yyyyrryyy', 'yyyyyryrr', 'yyyyyryry',
    'yyyyyryyr', 'yyyyyryyy', 'yyyyyyyyr', 'yyyyyyyyy'
]

# 转大写
patterns_upper = [p.upper() for p in patterns]

# 每 4 个一行输出
for i in range(0, len(patterns_upper), 4):
    group = patterns_upper[i:i+4]
    print(' '.join(group))
