# Doomsday Fuel
# =============

# Making fuel for the LAMBCHOP's reactor core is a tricky process because of the exotic matter involved.
# It starts as raw ore, then during processing, begins randomly changing between forms, eventually reaching
# a stable form. There may be multiple stable forms that a sample could ultimately reach,
# not all of which are useful as fuel.

# Commander Lambda has tasked you to help the scientists increase fuel creation efficiency by predicting
# the end state of a given ore sample. You have carefully studied the different structures that the ore
# can take and which transitions it undergoes. It appears that, while random, the probability of each
# structure transforming is fixed. That is, each time the ore is in 1 state, it has the same probabilities
# of entering the next state (which might be the same state).
# You have recorded the observed transitions in a matrix.
# The others in the lab have hypothesized more exotic forms that the ore can become,
# but you haven't seen all of them.

# Write a function solution(m) that takes an array of array of nonnegative ints representing how many times
# that state has gone to the next state and return an array of ints for each terminal state giving the exact
# probabilities of each terminal state, represented as the numerator for each state, then the denominator
# for all of them at the end and in simplest form. The matrix is at most 10 by 10. It is guaranteed that no
# matter which state the ore is in, there is a path from that state to a terminal state.
# That is, the processing will always eventually end in a stable state. The ore starts in state 0.
# The denominator will fit within a signed 32-bit integer during the calculation, as long as the fraction
# is simplified regularly.

# For example, consider the matrix m:
# [
#   [0,1,0,0,0,1],  # s0, the initial state, goes to s1 and s5 with equal probability
#   [4,0,0,3,2,0],  # s1 can become s0, s3, or s4, but with different probabilities
#   [0,0,0,0,0,0],  # s2 is terminal, and unreachable (never observed in practice)
#   [0,0,0,0,0,0],  # s3 is terminal
#   [0,0,0,0,0,0],  # s4 is terminal
#   [0,0,0,0,0,0],  # s5 is terminal
# ]
# So, we can consider different paths to terminal states, such as:
# s0 -> s1 -> s3
# s0 -> s1 -> s0 -> s1 -> s0 -> s1 -> s4
# s0 -> s1 -> s0 -> s5
# Tracing the probabilities of each, we find that
# s2 has probability 0
# s3 has probability 3/14
# s4 has probability 1/7
# s5 has probability 9/14
# So, putting that together, and making a common denominator, gives an answer in the form of
# [s2.numerator, s3.numerator, s4.numerator, s5.numerator, denominator] which is
# [0, 3, 2, 9, 14].

# Languages
# =========

# To provide a Java solution, edit Solution.java
# To provide a Python solution, edit solution.py

# Test cases
# ==========
# Your code should pass the following test cases.
# Note that it may also be run against hidden test cases not shown here.

# -- Java cases --
# Input:
# Solution.solution({{0, 2, 1, 0, 0}, {0, 0, 0, 3, 4}, {0, 0, 0, 0, 0}, {0, 0, 0, 0,0}, {0, 0, 0, 0, 0}})
# Output:
#     [7, 6, 8, 21]

# Input:
# Solution.solution({{0, 1, 0, 0, 0, 1}, {4, 0, 0, 3, 2, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}})
# Output:
#     [0, 3, 2, 9, 14]

# -- Python cases --
# Input:
# solution.solution([[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0,0], [0, 0, 0, 0, 0]])
# Output:
#     [7, 6, 8, 21]

# Input:
# solution.solution([[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
# Output:
#     [0, 3, 2, 9, 14]


from fractions import gcd
from fractions import Fraction


def maxtrixMultiplication(a, b):
    rows = len(a)
    cols = len(b[0])

    c = makeList(rows, cols)

    for row in range(rows):
        for col in range(cols):
            dotProduct = Fraction(0, 1)
            for i in range(len(a[0])):
                dotProduct += a[row][i]*b[i][col]
            c[row][col] = dotProduct
    return c


def multiplyMatrixRow(m, row, k):
    row_operator = makeIdentity(len(m))
    row_operator[row][row] = k
    return maxtrixMultiplication(row_operator, m)


def makeList(rows, cols):
    a = []
    for row in range(rows):
        a += [[0] * cols]
    return a


def makeIdentity(n):
    result = makeList(n, n)
    for i in range(n):
        result[i][i] = Fraction(1, 1)
    return result


def addMultipleSquareMatrix(m, srcRow, k, targetRow):
    row_operator = makeIdentity(len(m))
    row_operator[targetRow][srcRow] = k
    return maxtrixMultiplication(row_operator, m)


def invertMatrix(m):
    n = len(m)
    assert(len(m) == len(m[0]))
    inverse = makeIdentity(n)
    for col in range(n):
        diagonalRow = col
        assert(m[diagonalRow][col] != 0)
        k = Fraction(1, m[diagonalRow][col])
        m = multiplyMatrixRow(m, diagonalRow, k)
        inverse = multiplyMatrixRow(inverse, diagonalRow, k)
        srcRow = diagonalRow
        for targetRow in range(n):
            if srcRow != targetRow:
                k = -m[targetRow][col]
                m = addMultipleSquareMatrix(
                    m, srcRow, k, targetRow)
                inverse = addMultipleSquareMatrix(
                    inverse, srcRow, k, targetRow)
    return inverse


def transformMatrix(m):
    for rowIdx, row in enumerate(m):
        rowSum = sum(m[rowIdx])
        if rowSum == 0:
            m[rowIdx][rowIdx] = 1
        else:
            for colIdx, col in enumerate(row):
                m[rowIdx][colIdx] = Fraction(col, rowSum)


def getSubMatrix(m, rows, cols):
    newMatrix = []

    for row in rows:
        currRow = []
        for col in cols:
            currRow.append(m[row][col])
        newMatrix.append(currRow)
    return newMatrix


def get_q(m, nonTerminalStates):
    return getSubMatrix(m, nonTerminalStates, nonTerminalStates)


def get_r(m, nonTerminalStates, terminalStates):
    return getSubMatrix(m, nonTerminalStates, terminalStates)


def subtractMatrices(a, b):
    newMatrix = []
    for rowIdx, row in enumerate(a):
        column = []
        for colIdx, col in enumerate(row):
            column.append(a[rowIdx][colIdx] - b[rowIdx][colIdx])
        newMatrix.append(column)

    return newMatrix


def lcm(a, b):
    result = a * b / gcd(a, b)

    return result


def arrayLcm(args):
    arrLen = len(args)
    if arrLen <= 2:
        return lcm(*args)

    initial = lcm(args[0], args[1])
    i = 2
    while i < arrLen:
        initial = lcm(initial, args[i])
        i += 1
    return initial


def solution(m):
    terminalStates = []
    nonTerminalStates = []
    for index, row in enumerate(m):
        if sum(row) == 0:
            terminalStates.append(index)
        else:
            nonTerminalStates.append(index)

    if len(terminalStates) == 1:
        return [1, 1]

    transformMatrix(m)

    q = get_q(m, nonTerminalStates)
    r = get_r(m, nonTerminalStates, terminalStates)

    result = maxtrixMultiplication(invertMatrix(
        subtractMatrices(makeIdentity(len(q)), q)), r)

    denominator = arrayLcm([item.denominator for item in result[0]])

    result = [item.numerator * denominator /
              item.denominator for item in result[0]]

    result.append(denominator)

    return result


print(solution([[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [
      0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]) == [7, 6, 8, 21])

print(solution([[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [
      0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]) == [0, 3, 2, 9, 14])

# This ques took a lot from me. his was by far the most challenging problem for me. But this problem helped me
# to learn about markov chains and absorbing markov chains
