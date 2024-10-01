use std::{
    cmp::{max, Ordering},
    mem,
    ops::Not
};

struct AVLNode<T: PartialOrd> {
    value: T,
    height: usize,
    length: usize,
    left: Option<Box<AVLNode<T>>>,
    right: Option<Box<AVLNode<T>>>,
}

pub struct AVLTree<T: PartialOrd> {
    root: Option<Box<AVLNode<T>>>,
}

#[derive(Clone, Copy)]
enum Side {
    Left,
    Right,
}

impl<T: PartialOrd> AVLTree<T> {
    pub fn new() -> AVLTree<T> {
        AVLTree {
            root: None,
        }
    }

    pub fn insert(&mut self, value: T) -> usize {
        insert(&mut self.root, value)
    }
}

fn insert<T: PartialOrd>(tree: &mut Option<Box<AVLNode<T>>>, value: T) -> usize {
    if let Some(node) = tree {
        let count = match value.partial_cmp(&node.value) {
            Some(Ordering::Less) => insert(&mut node.left, value),
            _ => node.length(Side::Left) + 1 + insert(&mut node.right, value)
        };
        node.rebalance();
        count
    } else {
        *tree = Some(Box::new(AVLNode {
            value,
            height: 1,
            length: 1,
            left: None,
            right: None,
        }));
        0
    }
}

impl<T: PartialOrd> AVLNode<T> {
    fn child(&self, side: Side) -> &Option<Box<AVLNode<T>>> {
        match side {
            Side::Left => &self.left,
            Side::Right => &self.right,
        }
    }
    fn child_mut(&mut self, side: Side) -> &mut Option<Box<AVLNode<T>>> {
        match side {
            Side::Left => &mut self.left,
            Side::Right => &mut self.right,
        }
    }
    fn height(&self, side: Side) -> usize {
        self.child(side).as_ref().map_or(0, |n| n.height)
    }
    fn length(&self, side: Side) -> usize {
        self.child(side).as_ref().map_or(0, |n| n.length)
    }
    fn balance_factor(&self) -> i8 {
        let (left, right) = (self.height(Side::Left), self.height(Side::Right));
        if left < right { (right - left) as i8 } 
        else { -((left - right) as i8) }
    }
    fn update_height(&mut self) {
        self.height = 1 + max(self.height(Side::Left), self.height(Side::Right));
    }
    fn update_length(&mut self) {
        self.length = 1 + self.length(Side::Left) + self.length(Side::Right)
    }
    fn rotate(&mut self, side: Side) {
        let mut subtree = self.child_mut(!side).take().unwrap();
        *self.child_mut(!side) = subtree.child_mut(side).take();
        self.update_height();
        self.update_length();
        mem::swap(self, subtree.as_mut());
        *self.child_mut(side) = Some(subtree);
        self.update_height();
        self.update_length();
    }
    fn rebalance(&mut self) {
        self.update_height();
        self.update_length();
        let side = match self.balance_factor() {
            -2 => Side::Left,
            2 => Side::Right,
            _ => return,
        };
        let subtree = self.child_mut(side).as_mut().unwrap();
        if let (Side::Left, 1) | (Side::Right, -1) = (side, subtree.balance_factor()) {
            subtree.rotate(side);
        }
        self.rotate(!side);
    }
}

impl Not for Side {
    type Output = Side;

    fn not(self) -> Self::Output {
        match self {
            Side::Left => Side::Right,
            Side::Right => Side::Left,
        }
    }
}