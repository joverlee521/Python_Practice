'''
Next Prime Number
    A program that finds prime numbers until the user chooses to stop asking for the next one. 
'''
from math import floor
current_prime = 3

def get_next_prime(previous_prime):
    possible_prime = previous_prime + 2
    while True:
        isPrime = True
        for num in range(3, 3, 2):
            print("hi")
        
        for num in range(3, floor(possible_prime / 2), 2):
            if possible_prime % num == 0:
                isPrime = False
                possible_prime += 2
                break
        if isPrime:
            break
    return possible_prime


def main():
    global current_prime
    print("Let's look at some prime numbers!")
    print("The first prime number is 2.")
    while True:
        confirm = input("Do you want the next prime number(Y/N)?: ").upper()
        if confirm == "Y":
            print("The next prime number is {}".format(current_prime))
            current_prime = get_next_prime(current_prime)
            continue
        elif confirm == "N":
            break
        else: 
            print("Please answer with Y or N!")
            continue

if __name__ == "__main__":
    main()