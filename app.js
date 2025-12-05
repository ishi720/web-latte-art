"use strict;"

const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');

const N = 128;
const SIZE = (N + 2) * (N + 2);
const DT = 0.1;

// 流体シミュレーションの状態
let u = new Float32Array(SIZE);
let v = new Float32Array(SIZE);
let u_prev = new Float32Array(SIZE);
let v_prev = new Float32Array(SIZE);
let dens = new Float32Array(SIZE);
let dens_prev = new Float32Array(SIZE);

// コーヒーの背景色（ミルクが混ざっていない状態）
let coffee = new Float32Array(SIZE);

// マウス状態
let mouseDown = false;
let mouseX = 0, mouseY = 0;
let pmouseX = 0, pmouseY = 0;

// パラメータ
let viscosity = 0.005;
let diffusion = 0.000001;
let force = 5;
let densityAmount = 300;

function IX(i, j) {
    return i + (N + 2) * j;
}

function set_bnd(b, x) {
    const centerX = (N + 2) / 2;
    const centerY = (N + 2) / 2;
    const radius = N * 0.45;

    for (let i = 1; i <= N; i++) {
        x[IX(0, i)] = b === 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b === 1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = b === 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = b === 2 ? -x[IX(i, N)] : x[IX(i, N)];
    }

    // 円形の境界条件（カップの縁で跳ね返る）
    for (let j = 1; j <= N; j++) {
        for (let i = 1; i <= N; i++) {
            const dx = i - centerX;
            const dy = j - centerY;
            const dist = Math.sqrt(dx * dx + dy * dy);

            if (dist > radius) {
                // 境界の外側では速度を反転
                if (b === 1) { // u velocity
                    x[IX(i, j)] = 0;
                } else if (b === 2) { // v velocity
                    x[IX(i, j)] = 0;
                } else { // density
                    x[IX(i, j)] *= 0.5;
                }
            }
        }
    }

    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

function lin_solve(b, x, x0, a, c) {
    const cRecip = 1.0 / c;
    for (let k = 0; k < 20; k++) {
        for (let j = 1; j <= N; j++) {
            for (let i = 1; i <= N; i++) {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (
                    x[IX(i + 1, j)] + x[IX(i - 1, j)] +
                    x[IX(i, j + 1)] + x[IX(i, j - 1)]
                )) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}

function diffuse(b, x, x0, diff, dt) {
    const a = dt * diff * N * N;
    lin_solve(b, x, x0, a, 1 + 4 * a);
}

function advect(b, d, d0, u, v, dt) {
    const dt0 = dt * N;
    for (let j = 1; j <= N; j++) {
        for (let i = 1; i <= N; i++) {
            let x = i - dt0 * u[IX(i, j)];
            let y = j - dt0 * v[IX(i, j)];

            if (x < 0.5) x = 0.5;
            if (x > N + 0.5) x = N + 0.5;
            let i0 = Math.floor(x);
            let i1 = i0 + 1;

            if (y < 0.5) y = 0.5;
            if (y > N + 0.5) y = N + 0.5;
            let j0 = Math.floor(y);
            let j1 = j0 + 1;

            let s1 = x - i0;
            let s0 = 1 - s1;
            let t1 = y - j0;
            let t0 = 1 - t1;

            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                            s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(b, d);
}

function project(u, v, p, div) {
    for (let j = 1; j <= N; j++) {
        for (let i = 1; i <= N; i++) {
            div[IX(i, j)] = -0.5 * (
                u[IX(i + 1, j)] - u[IX(i - 1, j)] +
                v[IX(i, j + 1)] - v[IX(i, j - 1)]
            ) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 4);

    for (let j = 1; j <= N; j++) {
        for (let i = 1; i <= N; i++) {
            u[IX(i, j)] -= 0.5 * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
            v[IX(i, j)] -= 0.5 * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
        }
    }
    set_bnd(1, u);
    set_bnd(2, v);
}

function dens_step(x, x0, u, v, diff, dt) {
    for (let i = 0; i < SIZE; i++) {
        x[i] += dt * x0[i];
    }
    diffuse(0, x0, x, diff, dt);
    advect(0, x, x0, u, v, dt);
}

function vel_step(u, v, u0, v0, visc, dt) {
    for (let i = 0; i < SIZE; i++) {
        u[i] += dt * u0[i];
        v[i] += dt * v0[i];
    }
    diffuse(1, u0, u, visc, dt);
    diffuse(2, v0, v, visc, dt);
    project(u0, v0, u, v);
    advect(1, u, u0, u0, v0, dt);
    advect(2, v, v0, u0, v0, dt);
    project(u, v, u0, v0);
}

// マウスイベント
function getMousePos(e) {
    const rect = canvas.getBoundingClientRect();
    const scaleX = canvas.width / rect.width;
    const scaleY = canvas.height / rect.height;

    let clientX, clientY;
    if (e.touches) {
        clientX = e.touches[0].clientX;
        clientY = e.touches[0].clientY;
    } else {
        clientX = e.clientX;
        clientY = e.clientY;
    }

    return {
        x: (clientX - rect.left) * scaleX,
        y: (clientY - rect.top) * scaleY
    };
}

canvas.addEventListener('mousedown', (e) => {
    mouseDown = true;
    const pos = getMousePos(e);
    mouseX = pos.x;
    mouseY = pos.y;
    pmouseX = mouseX;
    pmouseY = mouseY;
});

canvas.addEventListener('mousemove', (e) => {
    pmouseX = mouseX;
    pmouseY = mouseY;
    const pos = getMousePos(e);
    mouseX = pos.x;
    mouseY = pos.y;

    if (mouseDown) {
        addMilk();
    }
});

canvas.addEventListener('mouseup', () => {
    mouseDown = false;
});

canvas.addEventListener('mouseleave', () => {
    mouseDown = false;
});

canvas.addEventListener('touchstart', (e) => {
    e.preventDefault();
    mouseDown = true;
    const pos = getMousePos(e);
    mouseX = pos.x;
    mouseY = pos.y;
    pmouseX = mouseX;
    pmouseY = mouseY;
});

canvas.addEventListener('touchmove', (e) => {
    e.preventDefault();
    pmouseX = mouseX;
    pmouseY = mouseY;
    const pos = getMousePos(e);
    mouseX = pos.x;
    mouseY = pos.y;

    if (mouseDown) {
        addMilk();
    }
});

canvas.addEventListener('touchend', () => {
    mouseDown = false;
});

// ミルクを注ぐ
function addMilk() {
    const i = Math.floor((mouseX / canvas.width) * N) + 1;
    const j = Math.floor((mouseY / canvas.height) * N) + 1;

    if (i >= 1 && i <= N && j >= 1 && j <= N) {
        const dx = mouseX - pmouseX;
        const dy = mouseY - pmouseY;

        // 速度を追加
        const radius = 3;
        for (let di = -radius; di <= radius; di++) {
            for (let dj = -radius; dj <= radius; dj++) {
                const ni = i + di;
                const nj = j + dj;
                if (ni >= 1 && ni <= N && nj >= 1 && nj <= N) {
                    const dist = Math.sqrt(di * di + dj * dj);
                    if (dist <= radius) {
                        const factor = (1 - dist / radius);
                        u_prev[IX(ni, nj)] += dx * force * factor;
                        v_prev[IX(ni, nj)] += dy * force * factor;
                        dens_prev[IX(ni, nj)] += densityAmount * factor;
                    }
                }
            }
        }
    }
}

// カフェアート風の描画
function draw() {
    const imageData = ctx.createImageData(canvas.width, canvas.height);
    const data = imageData.data;

    const cellWidth = canvas.width / N;
    const cellHeight = canvas.height / N;

    const centerX = N / 2;
    const centerY = N / 2;
    const cupRadius = N * 0.45;

    for (let j = 1; j <= N; j++) {
        for (let i = 1; i <= N; i++) {
            const dx = i - centerX;
            const dy = j - centerY;
            const distFromCenter = Math.sqrt(dx * dx + dy * dy);

            // カップの外側は描画しない
            if (distFromCenter > cupRadius) {
                continue;
            }

            const d = dens[IX(i, j)];
            const c = coffee[IX(i, j)];

            // ミルク（白）とコーヒー（茶色）をブレンド
            const milkAmount = Math.min(1, d / 100);

            // コーヒー色: 濃い茶色
            const coffeeR = 40 + c * 0.3;
            const coffeeG = 25 + c * 0.2;
            const coffeeB = 15 + c * 0.1;

            // ミルク色: クリーム色
            const milkR = 245;
            const milkG = 235;
            const milkB = 220;

            // ブレンド
            let r = coffeeR * (1 - milkAmount) + milkR * milkAmount;
            let g = coffeeG * (1 - milkAmount) + milkG * milkAmount;
            let b = coffeeB * (1 - milkAmount) + milkB * milkAmount;

            // カップの縁に向かって暗くする（ビネット効果）
            const vignette = 1 - (distFromCenter / cupRadius) * 0.3;
            r *= vignette;
            g *= vignette;
            b *= vignette;

            r = Math.min(255, Math.max(0, r));
            g = Math.min(255, Math.max(0, g));
            b = Math.min(255, Math.max(0, b));

            // ピクセルを塗りつぶす
            const x = (i - 1) * cellWidth;
            const y = (j - 1) * cellHeight;

            for (let py = 0; py < cellHeight; py++) {
                for (let px = 0; px < cellWidth; px++) {
                    const pixelX = Math.floor(x + px);
                    const pixelY = Math.floor(y + py);

                    // カップの円形マスク
                    const pdx = pixelX - canvas.width / 2;
                    const pdy = pixelY - canvas.height / 2;
                    const pixelDist = Math.sqrt(pdx * pdx + pdy * pdy);

                    if (pixelX < canvas.width && pixelY < canvas.height &&
                        pixelDist <= canvas.width / 2) {
                        const idx = (pixelY * canvas.width + pixelX) * 4;
                        data[idx] = r;
                        data[idx + 1] = g;
                        data[idx + 2] = b;
                        data[idx + 3] = 255;
                    }
                }
            }
        }
    }

    ctx.putImageData(imageData, 0, 0);
}

// メインループ
function simulate() {
    vel_step(u, v, u_prev, v_prev, viscosity, DT);
    dens_step(dens, dens_prev, u, v, diffusion, DT);

    u_prev.fill(0);
    v_prev.fill(0);
    dens_prev.fill(0);

    // 速度を減衰させて、数秒で拡散が止まるようにする
    for (let i = 0; i < SIZE; i++) {
        u[i] *= 0.97;
        v[i] *= 0.97;
    }

    // 密度の減衰を非常に少なくする（ミルクのパターンが残る）
    for (let i = 0; i < SIZE; i++) {
        dens[i] *= 0.9995;
    }

    draw();
    requestAnimationFrame(simulate);
}

// コントロール
document.getElementById('viscosity').addEventListener('input', (e) => {
    viscosity = e.target.value / 100000;
    document.getElementById('viscosity-value').textContent = viscosity.toFixed(6);
});

document.getElementById('diffusion').addEventListener('input', (e) => {
    diffusion = e.target.value / 1000000;
    document.getElementById('diffusion-value').textContent = diffusion.toFixed(6);
});

document.getElementById('force').addEventListener('input', (e) => {
    force = parseInt(e.target.value);
    document.getElementById('force-value').textContent = force;
});

document.getElementById('density').addEventListener('input', (e) => {
    densityAmount = parseInt(e.target.value);
    document.getElementById('density-value').textContent = densityAmount;
});

// プリセット機能
function resetCoffee() {
    u.fill(0);
    v.fill(0);
    u_prev.fill(0);
    v_prev.fill(0);
    dens.fill(0);
    dens_prev.fill(0);

    // コーヒーの背景をランダムに
    for (let i = 0; i < SIZE; i++) {
        coffee[i] = Math.random() * 20;
    }
}

// 初期化
resetCoffee();
simulate();