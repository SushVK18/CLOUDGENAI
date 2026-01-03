const jwt = require('jsonwebtoken');
const jwksClient = require('jwks-rsa');
const audience = process.env.AUDIENCE;
const algorithms = process.env.ALGORITHMS || 'RS256';
const issuersUrls = process.env.ISSUERSURLS || 'https://sts.windows.net/fcb2b37b-5da0-466b-9b83-0014b67a7c78/';
const issuers = issuersUrls.split(',');

const verificationOptions = {
    algorithms,
    ...(audience ? { audience } : {}),
    issuer: issuers
};

function generatePolicy(principalId, effect, resource) {
    const authResponse = {
        principalId,
        policyDocument: {
            Version: '2012-10-17',
            Statement: [{
                Action: 'execute-api:Invoke',
                Effect: effect,
                Resource: resource
            }]
        }
    };
    console.log("Generated policy response:", JSON.stringify(authResponse));
    return authResponse;
}

function getJwksUri(token) {
    const payload = token.split('.')[1];
    let decodedPayload;
    try {
        decodedPayload = JSON.parse(Buffer.from(payload, 'base64').toString());
    } catch (error) {
        console.error("Error parsing JWT payload:", error);
        throw new Error("Unauthorized");
    }
    return decodedPayload.ver === '2.0' ? 'https://login.microsoftonline.com/common/discovery/v2.0/keys' : 'https://login.microsoftonline.com/common/discovery/keys';
}

function getSigningKey(header, callback, keyClient) {
    keyClient.getSigningKey(header.kid, function (err, key) {
        if (err) {
            callback(err, null);
        } else {
            const signingKey = key.publicKey || key.rsaPublicKey;
            callback(null, signingKey);
        }
    });
};


function validateToken(token, methodArn, callback) {
    const keyClient = jwksClient({
        cache: true,
        cacheMaxAge: 600000,
        strictSsl: true,
        jwksUri: getJwksUri(token)
    });
    jwt.verify(token, (header, callback) => getSigningKey(header, callback, keyClient), verificationOptions, function (error, decoded) {
        console.log("Decoded:", decoded);
        console.log("Expected issuers:", issuers);
        console.log("Expected audience:", audience);
        if (error) {
            console.error("JWT verification error:", error);
            callback("Unauthorized");
        } else {
            if (!issuers.includes(decoded.iss)) {
                console.error("Invalid issuer:", decoded.iss);
                callback("Unauthorized");
            } else {
                callback(null, generatePolicy(decoded.sub, "Allow", methodArn));
            }
        }
    });
}

exports.handler = (event, context, callback) => {
    const authHeader = event.headers.Authorization || event.headers.authorization;
    if (!authHeader || !authHeader.toLowerCase().startsWith("bearer ")) {
      return callback(null, { isAuthorized: false });
    }
  
    const token = authHeader.slice(7);
    const keyClient = jwksClient({ jwksUri: getJwksUri(token) });
  
    jwt.verify(token, (header, cb) => {
      keyClient.getSigningKey(header.kid, (err, key) => {
        if (err) cb(err); else cb(null, key.getPublicKey());
      });
    }, { audience, algorithms: [algorithms], issuer: issuers }, (err, decoded) => {
      if (err) {
        console.log("JWT error:", err);
        return callback(null, { isAuthorized: false });
      }
  
      return callback(null, {
        isAuthorized: true,
        context: {
          sub: decoded.sub,
          aud: decoded.aud,
          iss: decoded.iss
        }
      });
    });
  };
