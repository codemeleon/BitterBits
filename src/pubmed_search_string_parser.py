#!/usr/bin/env python

"""Converstion of expresion to tree."""



# Function to check parentheses
def check_balance(myStr):
    """Check the string is Balanced or not"""
    open_list = ["[","{","("]
    close_list = ["]","}",")"]

    stack = []
    for i in myStr:
        if i in open_list:
            stack.append(i)
        elif i in close_list:
            pos = close_list.index(i)
            if ((len(stack) > 0) and
                (open_list[pos] == stack[len(stack)-1])):
                stack.pop()
            else:
                return "Unbalanced"
    if len(stack) == 0:
        return 1 # Balanced
    else:
        return 0 # Unbalanced




class Node:

    """Docstring for Node. """

    def __init__(self, data):
        """TODO: to be defined.

        :data: TODO

        """

        self.left = None
        self.right = None
        self.data = data
        self.boolean = 0  # -1:or, 1: & and 0: not

    def insert(self, data):
        """TODO: Docstring for insert.

        :data: TODO
        :returns: TODO

        """
        if self.data
        pass








def string2tree(
        qstring="((Anmol[Author]) OR (RNA Editind)) AND \
                (Test Subject[Author])"):
    """This converts string to binary three

    :qstring: TODO
    :returns: TODO
    TODO: Check if variable are allowed for the search, if not put in all or
    throw error
    ~ : Not
    | : or
    & : and

    """
    pass


def tree2dquery(tree):
    """Converts tree to django query and return the query string

    :tree: TODO
    :returns: TODO

    """
    pass
