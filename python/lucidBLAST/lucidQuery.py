__author__ = 'Hunter Cameron'

from lucidParent import Parent

# TODO Recondsider the essential nature of this class and its imperitive single statement...
class Query(Parent):

    def __init__(self, args):
        Parent.__init__(self, args["name"], args["length"])

