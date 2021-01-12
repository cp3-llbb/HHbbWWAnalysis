import tensorflow as tf

@tf.function
def grouped_cross_entropy_t(
    labels,
    predictions,
    sample_weight=None,
    group_ids=None,
    focal_gamma=None,
    class_weight=None,
    epsilon=1e-7,
):
    assert group_ids is not None
    # get true-negative component
    predictions = tf.clip_by_value(predictions, epsilon, 1 - epsilon)
    tn = labels * tf.math.log(predictions)
    # focal loss?
    if focal_gamma is not None:
        tn *= (1 - predictions) ** focal_gamma
    # convert into loss
    losses = -tn
    # apply class weights
    if class_weight is not None:
        losses *= class_weight
    # apply sample weights
    if sample_weight is not None:
        losses *= sample_weight[:, tf.newaxis]
    # create grouped labels and predictions
    labels_grouped = []
    for _, ids in group_ids:
        labels_grouped.append(
            tf.reduce_sum(tf.gather(labels, ids, axis=-1), axis=-1, keepdims=True)
        )
    labels_grouped = (
        tf.concat(labels_grouped, axis=-1) if len(labels_grouped) > 1 else labels_grouped[0]
    )
    predictions_grouped = []
    for _, ids in group_ids:
        predictions_grouped.append(
            tf.reduce_sum(tf.gather(predictions, ids, axis=-1), axis=-1, keepdims=True)
        )
    predictions_grouped = (
        tf.concat(predictions_grouped, axis=-1)
        if len(predictions_grouped) > 1
        else predictions_grouped[0]
    )

    predictions_grouped = tf.clip_by_value(predictions_grouped, epsilon, 1 - epsilon)
    # grouped true-negative component
    tn_grouped = labels_grouped * tf.math.log(predictions_grouped)
    # focal loss?
    if focal_gamma is not None:
        tn_grouped *= (1 - predictions_grouped) ** focal_gamma
    # convert into loss and apply group weights
    group_weights = tf.constant([w for w, _ in group_ids], tf.float32)
    losses_grouped = -tn_grouped * group_weights
    # apply sample weights
    if sample_weight is not None:
        losses_grouped *= sample_weight[:, tf.newaxis]
    # combine losses
    loss = tf.reduce_mean(
        0.5 * (tf.reduce_sum(losses, axis=-1) + tf.reduce_sum(losses_grouped, axis=-1))
    )
    return loss


# Custom Loss Functions


class GroupedXEnt(tf.keras.losses.Loss):
    def __init__(
        self,
        group_ids=None,
        focal_gamma=None,
        class_weight=None,
        epsilon=1e-7,
        *args,
        **kwargs,
    ):
        super(GroupedXEnt, self).__init__(*args, **kwargs)
        self.group_ids = group_ids
        self.focal_gamma = focal_gamma
        self.class_weight = class_weight
        self.epsilon = epsilon

    def call(self, y_true, y_pred, sample_weight=None):
        return grouped_cross_entropy_t(
            y_true,
            y_pred,
            sample_weight=sample_weight,
            group_ids=self.group_ids,
            focal_gamma=self.focal_gamma,
            class_weight=self.class_weight,
            epsilon=self.epsilon,
        )


