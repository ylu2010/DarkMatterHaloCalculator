from wtforms import Form, FloatField, validators

class InputForm(Form):
    Omega_L = FloatField(
        label='\( \Omega_{\Lambda,0} \) ', description=' ', default=0.7,
        validators=[validators.InputRequired()])
    Omega_M = FloatField(
        label='\( \Omega_{M,0} \) ', description=' ', default=0.3,
        validators=[validators.InputRequired()])
    Hubble = FloatField(
        label='Hubble', description=' \({km}\,{s}^{-1}\,{Mpc}^{-1} \)', default=70,
        validators=[validators.InputRequired()])
    Mvir = FloatField(
        label='Mass', description=' \( M_{\odot} \)', default=1e12,
        validators=[validators.InputRequired()])
    Redshift = FloatField(
        label='Redshift', description=' ', default=0.0,
        validators=[validators.InputRequired()])

    #r = FloatField(validators=[validators.InputRequired()])
