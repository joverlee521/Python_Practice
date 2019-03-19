'''

Prime Factorization
    Have the user enter a number and find all Prime Factors (if there are any) and display them.

'''
from math import floor

def generate_primes(n):
    '''

    Generates a prime number sequence that is less than the given integer
    
    INPUT: n(int): the upper limit of the prime number sequence

    '''
    for num in range(2, n):
        is_prime = True
        # immediately yield 2 as prime
        if num == 2:
            yield num
        # check if num is an even number, which is definitely not prime
        elif num != 2 and num % 2 == 0:
            is_prime = False
        else:
            # Check if num is divisible by all odd numbers less than its square root to find out if it's a prime number
            for i in range(3, floor(num ** 0.5), 2):
                if num % i == 0:
                    is_prime = False
        if is_prime:
            yield num

def find_prime_factors(number):
    '''

    Generates a prime factor sequence for the given integer

    INPUT: number(int): number for findind prime factors 

    '''
    for prime_number in generate_primes(number):
        # Find prime numbers that are factors of given number
        while number % prime_number == 0:
            # Continue dividing number by current prime number until no longer divisible
            number = number / prime_number
            yield prime_number
        # If the number is equal to 1, then we have found all of the prime factors
        if number == 1:
            break

def main():
    '''

    Asks user to enter a number to find all of its prime factors
    and prints out a list of the prime numbers.

    If there are no prime factors, informs user their number is a prime number

    '''
    result = []
    while True:
        try: 
            number = int(input("Enter a number to find all prime factors: "))
        except: 
            print("Please enter a positive integer!")
            continue
        else:
            if number < 0:
                print("Please enter a positive integer!")
                continue
            for num in find_prime_factors(number):
                result.append(num)
            break
    if len(result) > 0:
        return print(result)
    else:
        return print(f"{number} is a prime number!")

if __name__ == "__main__":
    main()