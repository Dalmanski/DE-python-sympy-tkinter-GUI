def parse_user_input(user_input):
    user_data = user_input.split(',')

    p = None
    rp = None
    rt = None
    t = None

    for entry in user_data:
        variable_name, value = entry.split('=')

        if variable_name == 'p':
            p = value
        elif variable_name == 'rp':
            rp = value
        elif variable_name == 'rt':
            rt = value
        elif variable_name == 't':
            t = value
        else:
            print(f"variable {variable_name} is not define.")

    return p, rp, rt, t

user_input = "p=5000, rp=15, rts=10, t=30"
user_input = user_input.replace(' ','')
p, rp, rt, t = parse_user_input(user_input)

print(f"p: {p}")
print(f"rp: {rp}")
print(f"rt: {rt}")
print(f"t: {t}")
