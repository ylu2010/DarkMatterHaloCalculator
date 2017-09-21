from flask import Flask, render_template, request
from compute import compute
from model import InputForm
import parameter as par

#global Omega_M, Omega_L, Hubble


import cosmology_lib as Cosmo
reload(Cosmo)


app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm(request.form)
    if request.method == 'POST' and form.validate():
        par.Omega_M = form.Omega_M.data
        par.Omega_L = form.Omega_L.data
        par.Hubble = form.Hubble.data
        reload(Cosmo)
        m = form.Mvir.data
        z = form.Redshift.data
        rvir, vvir, tvir, result = compute(m, z)
        #rvir = m
    else:
        rvir = None
        vvir = None
        tvir = None
        result = None

    return render_template("view1.html", form=form, rvir=rvir, vvir=vvir, tvir=tvir, result=result)


if __name__ == '__main__':
    app.run(debug=True)
