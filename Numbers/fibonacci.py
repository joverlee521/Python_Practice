'''

Fibonacci Sequence
    Enter a number and have the program generate the Fibonacci sequence to that number or to the Nth number.

'''

def fibonacci(n):
    '''
    Generates the Fibonacci Sequence with a size of n

    INPUT: n(int): the size of the Fibonacci Sequence
    '''
    first_num = 0
    second_num = 1
    for _ in range(n):
        yield first_num
        first_num, second_num = second_num, first_num + second_num

def main():
    '''
    Asks user to enter the size of the Fibonacci Sequence they want
    and prints out a list of the Fibonacci Sequence

    '''
    result = []
    while True: 
        try:
            fib_num = int(input("How many numbers of the Fibonacci Series do you want?: "))
        except:
            print("Please enter a positive integer!")
            continue
        else:
            for number in fibonacci(fib_num):
                result.append(number)
            break
    return print(result)

if __name__== "__main__":
    main()