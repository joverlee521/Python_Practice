class Node:
    def __init__(self, val = None, next_node = None):
        self.val = val
        self.next_node = next_node

    def __str__(self):
        return "My value is {}".format(self.val)

class Linked_List:
    def __init__(self):
        self.head = None
        self.size = 0

    def get(self, index):
        if index >= self.size or index < 0:
            return print("Invalid Index!")
        current = self.head
        for _ in range(index):
            current = current.next_node
        return current.val
    
    def addHead(self, val):
        new_next = None
        if self.head is not None:
            new_next = self.head
        self.head = Node(val, new_next)
        self.size += 1
    
    def addAtIndex(self, index, val):
        if index < 0 or index > self.size:
            return print("Invalid Index!")
        if index == 0:
            return self.addHead(val)
        if index == self.size:
            return self.addTail(val)
        current = self.head
        for _ in range(1, index):
            current = current.next_node
        new_next = current.next_node
        current.next_node = Node(val, new_next)
        self.size += 1

    def addTail(self, val):
        if self.head is None:
            self.head = Node(val)
        else:
            current = self.head
            while current.next_node is not None:
                current = current.next_node
            current.next_node = Node(val)
        self.size += 1

    def deleteAtIndex(self, index):
        if index < 0 or index >= self.size:
            return print("Invalid Index!")
        if index == 0:
            to_delete = self.head
            self.head = self.head.next_node
            del to_delete
        else:
            current = self.head
            node_before = None
            for i in range(index):
                if i == (index - 1):
                    node_before = current
                current = current.next_node
            node_before.next_node = current.next_node
            del current
        self.size -= 1
