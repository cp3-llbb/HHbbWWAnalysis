import os
import sys
import json
import argparse
import talos
import numpy as np
import tensorflow as tf
from lbn import LBNLayer
import Operations

class WhereEquals(tf.keras.layers.Layer):
    def __init__(self, value=0):
        super(WhereEquals, self).__init__()
        self.value = value
    def call(self, inp):
        return tf.where(inp[:, 0] == self.value)

def import_tf():
    import tensorflow as tf
    # keep a reference to the v1 API as long as v2 provides compatibility
    tf1 = None
    tf_version = tf.__version__.split(".", 2)
    if tf_version[0] == "1":
        tf1 = tf
    elif getattr(tf, "compat", None) is not None and getattr(tf.compat, "v1", None) is not None:
        tf1 = tf.compat.v1

    return tf, tf1, tf_version

class KerasToTensorflowModel:
    def __init__(self,name,outdir,paths_zip,paths_json,paths_h5):
        self.name       = name
        self.outdir     = outdir
        if paths_zip is None and paths_json is None and paths_h5 is None:
            raise RuntimeError('Either zip paths or json+h5 paths need to be provided')
        if paths_zip is not None and (paths_json is not None or paths_h5 is not None):
            raise RuntimeError('Zip paths cannot be used with either json or h5 paths')
        if paths_zip is None and not (paths_json is not None and paths_h5 is not None):
            raise RuntimeError('Json and h5 paths need to be used together')

        if paths_zip is not None:
            models = [self.load_model_talos(path_zip) for path_zip in paths_zip] 
        else:
            models = [self.load_model_keras(path_json,path_h5) for path_json,path_h5 in zip(paths_json,paths_h5)]

        if len(models) == 1:
            self.save_model_to_protobuf(models[0])
            self.makeInputList(models[0])
        elif len(models) > 1 :
            print ("Will stitch models :")
            for i in range(len(models)):
                if paths_zip is not None:
                    path = paths_zip[i]
                else:
                    path = paths_json[i]
                print ("... model {} applied on events with eventnr % {} == {}".format(path,len(models),i))
            stitch_model = self.stich_models(models)
            if self.test_stich(stitch_model,models):
                self.save_model_to_protobuf(stitch_model)
                self.makeInputList(stitch_model)
            else:
                raise RuntimeError("Stitching test failed")
        else:
            raise RuntimeError("No model has been loaded")
            

    @staticmethod
    def load_model_keras(path_to_json,path_to_h5):
        with open(path_to_json,"r") as f:
            loaded_model_json = f.read()
        custom_objects = {}

        if 'LBNLayer' in loaded_model_json:
            custom_objects.update({'LBNLayer':LBNLayer})
        if 'CategoryEncoding' in loaded_model_json:
            custom_objects.update({name:getattr(Operations,name) for name in dir(Operations) if name.startswith('op')})

        model = tf.keras.models.model_from_json(loaded_model_json,custom_objects=custom_objects)
        model.load_weights(path_to_h5)
        return model
            
    @staticmethod
    def load_model_talos(path_to_zip):
        custom_objects = {'LBNLayer':LBNLayer}
        custom_objects.update({name:getattr(Operations,name) for name in dir(Operations) if name.startswith('op')})
        model = talos.Restore(path_to_zip,custom_objects=custom_objects).model
        return model

    @staticmethod
    def stich_models(models):
        inputs = [tf.keras.Input(inp.shape[1:], name=inp.name) for inp in models[0].inputs]
        eventnr = tf.keras.Input((1,), dtype=tf.int64, name="eventnr")
        const = len(models)
        model_idx = tf.keras.layers.Lambda(lambda x: x % const)(eventnr)

        output = tf.zeros(
            shape=tf.concat((tf.shape(inputs[0])[:1], tf.shape(models[0](inputs))[1:]), axis=0)
        )
        for i, model in enumerate(models):
            #print ("Stitching model {}".format(i),end=' ... ')
            model._name += str(i)
            idx = WhereEquals(value=i)(model_idx)
            inp_gathered = [tf.gather_nd(a, idx) for a in inputs]
            out = model(inp_gathered)
            output = tf.tensor_scatter_nd_update(output, idx, out)
            print ('done')

        stitched_model = tf.keras.Model(inputs + [eventnr], output)

        return stitched_model

    @staticmethod
    def test_stich(stitched_model,models):
        batch = 100
        success = True
        for idx,model in enumerate(models):
            print ('Testing model {}'.format(idx),end=" ... ")
            input_list = [np.random.rand(batch,*tuple(inp.shape[1:])) for inp in models[0].inputs]
            output_true = model.predict(input_list)
            input_list.append(np.ones((batch,1))*idx)
            output_stitched = stitched_model.predict(input_list)
            if not (output_true == output_stitched).all():
                print("model not producing same output as stitched model")
                success = False    
            else:
                print ("success")
                
        return success


    def save_model_to_protobuf(self,model):
 #       model_path = os.path.join(self.outdir,self.name)
 #       model.save(model_path)
 #       print ("Saved model as {}".format(model_path))

        # convert keras models and polymorphic functions to concrete functions

        tf, tf1, tf_version = import_tf()

        assert tf_version[0] != "1"

        from tensorflow.python.keras.saving import saving_utils
        from tensorflow.python.eager.def_function import Function
        from tensorflow.python.eager.function import ConcreteFunction

        assert isinstance(model, tf.keras.Model)
        learning_phase_orig = tf.keras.backend.learning_phase()
        tf.keras.backend.set_learning_phase(False)
        model_func = saving_utils.trace_model_call(model)
        if model_func.function_spec.arg_names and not model_func.input_signature:
            raise ValueError("when model is a keras model callable accepting arguments, its "
                             "input signature must be frozen by building the model")
        model = model_func.get_concrete_function()
        tf.keras.backend.set_learning_phase(learning_phase_orig)

        # convert variables to constants
        assert isinstance(model, ConcreteFunction)

        from tensorflow.python.framework import convert_to_constants

        model = convert_to_constants.convert_variables_to_constants_v2(model)
        graph = model.graph

        tf.io.write_graph(graph, self.outdir, self.name+'.pb',as_text=False)
        tf.io.write_graph(graph, self.outdir, self.name+'.pb.txt',as_text=True)
        print ("Saved model as {}".format(os.path.join(self.outdir,self.name+'.pb')))

    def makeInputList(self,model):
        txtfile = os.path.join(self.outdir,self.name+'_inputs.txt')
        with open(txtfile,"w") as f:
            for layer in model.layers:
                if 'InputLayer' in  str(type(layer)):
                    f.write("{:50} -> {}\n".format(layer.name,str(layer.input_shape)))

        print ("Saved input txt file as {}".format(txtfile))
                    
         
        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser('tf.keras model converter to protobuf format\nNote that in case of stitched model (when doing cross-validation) several models can be provided but no check can be done on the correct ordering')
    parser.add_argument('--zip',required=False, type=str, nargs='+', default=None,
                        help='The zip model(s) file you wish to convert to .pb')
    parser.add_argument('--json',required=False, type=str, nargs='+', default=None,
                        help='The json model(s) file you wish to convert to .pb')
    parser.add_argument('--h5',required=False, type=str,  nargs='+', default=None,
                        help='The h5 model(s) model weights file you wish to convert to .pb')
    parser.add_argument('--outdir','-o', dest='outdir', required=False, default='TFModels/', 
                        help='The directory to place the output files - default("TFModels/")')
    parser.add_argument('--name', required=False, default='output_graph', 
                        help='The name of the resulting output - default("output_graph") (no .pb)')
    parser.add_argument('--test', required=False, default=None,
                        help='Test pb file')
    args = parser.parse_args()

    if args.test:
        tf, tf1, tf_version = import_tf()

        graph = tf.Graph()
        with graph.as_default():
            graph_def = graph.as_graph_def()

            # as text #
            if args.test.endswith("txt"):
                print ('using txt file')
                from google.protobuf import text_format
                with open(args.test, "rb") as f:
                    text_format.Merge(f.read(), graph_def)
            else:
                print ('using pb file')
                from google.protobuf import text_format
                with tf.io.gfile.GFile(args.test, "rb") as f:
                    graph_def.ParseFromString(f.read())
                
            tf.import_graph_def(graph_def, name="")

        sys.exit()


    instance = KerasToTensorflowModel(name          = args.name, 
                                      outdir        = args.outdir,
                                      paths_zip     = args.zip,
                                      paths_json    = args.json,
                                      paths_h5      = args.h5)
