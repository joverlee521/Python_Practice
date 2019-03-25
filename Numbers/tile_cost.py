'''

Find Cost of Tile to Cover W x H Floor
    Calculate the total cost of tile it would take to cover a floor plan of width and height, using a cost entered by the user.

'''
import decimal

def calc_cost(list_of_variables):
    result = 1
    for value in list_of_variables:
        result *= value
    return result

def main():
    outputs = {0: "width of the floor", 1:"length of the floor", 2:"cost per square meter"}
    variables = ["width", "length", "cost"]
    for key in outputs:
        try: 
            variables[key] = float(input(f"What is the {outputs[key]}?: "))
        except:
            print("Numbers only!")
            continue
        else: 
            if variables[key] <= 0: 
                print("Positive numbers only!")
                continue
    return print("The total cost is ${}".format(decimal.Decimal(calc_cost(variables)).quantize(decimal.Decimal(".01"))))
    

if __name__ == "__main__":
    main()