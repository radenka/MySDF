class Man:
    "Ideal partner evaluation. 100% reliable."

    def __init__(self, salary, your_salary): # (self, age, car, salary, debts)
        self.salary = salary
        self.your_salary = your_salary

    def count_money(self):
        if self.salary <= self.your_salary:
            print('No money. Dump him.')
        elif self.salary <= 1.5*self.your_salary:
            print('Not very suitable for you. Little money, no perspective.')
        elif self.salary <= 3*self.your_salary:
            print("Your target is capable of buying you new clothes every week.\nVery suitable for you.")
        else:
            print('Marry him. NOW.')

Marek = Man(13400, 1600)
Marek.count_money()

# Marek.salary = 10000  # changing values in object instance with 'instance.method = new_value'
# Marek.your_salary = 10000
# Marek.count_money()

