from flask import Flask
from flask_cors import CORS   # ✨ NEW
from routes.predict import predict_bp

app = Flask(__name__)

# ✅ Enable CORS for all routes (frontend -> backend communication)
CORS(app)

# ✅ Register your routes
app.register_blueprint(predict_bp, url_prefix='/api')

@app.route('/')
def home():
    return {'message': 'ChemStruct AI Backend is running!'}

if __name__ == '__main__':
    # run on default port 5000
    app.run(debug=True)
