# driver class. We'll expose the functions of this class
# after it's been instantiated

class ExampleDriver():

    def __init__(self, info):
        self.info = info

    def do_something(self):
        print("I'm alive!")

    def close(self):
        print("bye!")