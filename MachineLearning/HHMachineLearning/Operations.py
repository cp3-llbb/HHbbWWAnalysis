import tensorflow as tf


class OperationBase(tf.keras.layers.Layer):
    def __init__(self,**kwargs):
        super(OperationBase,self).__init__(kwargs)
    def call(self,x):
        return self.operation(x)
    def operation(self,x):
        raise NotImplementedError
    def compute_output_shape(self,input_shape):
        return input_shape
    @property
    def onehot_dim(self):
        raise NotImplementedError

class op_pdgid(OperationBase):
    def operation(self,x):
        return tf.cast(tf.abs(x) > 11 , "int32")
    @property
    def onehot_dim(self):
        return 2

class op_charge(OperationBase):
    def operation(self,x):
        return tf.cast(tf.abs(x) > 0, "int32")
    @property
    def onehot_dim(self):
        return 2

class op_era(OperationBase):
    def operation(self,x):
        return tf.cast(x - 2016 , "int32")
    @property
    def onehot_dim(self):
        return 3

class op_resolved_JPAcat(OperationBase):
    def operation(self,x):
        return tf.cast(x - 1, "int32")
    @property
    def onehot_dim(self):
        return 7

class op_boosted_JPAcat(OperationBase):
    def operation(self,x):
        return tf.cast(x - 1, "int32")
    @property
    def onehot_dim(self):
        return 3
