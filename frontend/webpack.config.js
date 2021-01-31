const path = require('path');
const fs = require('fs');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const {CleanWebpackPlugin} = require('clean-webpack-plugin');
const ZipWebpackPlugin = require('zip-webpack-plugin');
const CopyWebpackPlugin = require('copy-webpack-plugin');
const marked = require("marked");

let MODE = 'production';

if (process.argv.length === 3) {
    let mode_arg = process.argv[2].split('=')[0];
    let mode_value = process.argv[2].split('=')[1];
    if (mode_arg === '--mode' && mode_value === 'development') {
        MODE = 'development';
    }
}
DIST_PATH = path.resolve(__dirname, 'dist');

// Generate about HTML from markdown file
ABOUT_FILE = 'doc/about.md';
aboutMdHtml = marked(fs.readFileSync(ABOUT_FILE, 'utf8'));

aboutMdHtml = aboutMdHtml.replace('<table>', '<table class="table">');

module.exports = {
    mode: MODE,
    entry: {
        lrr_db: "./apps/lrr_db/lrr_db.js",
        find_lrr: "./apps/find_lrr/find_lrr.js",
    },
    output: {
        path: DIST_PATH,
        publicPath: "./",
        filename: MODE === 'production' ? 'js/[name].[contenthash].js' : 'js/[name].js',
    },
    externals: {
        jquery: 'jQuery',
        vue: 'Vue',
        'vue-multiselect': 'Multiselect'
    },
    module: {
        rules: [
            {
                test: /\.ejs$/,
                loader: 'ejs-loader'
            },
            {
                test: /\.m?js$/,
                exclude: /node_modules/,        //排除node_modules文件夹下的js
                use: {
                    loader: 'babel-loader',     //使用babel-loader处理找到的js文件
                    options: {                  ////采用babel-loader的"es2015"规则将找的js为浏览器可识别的js
                        presets: ['@babel/preset-env'],
                        plugins: ['@babel/plugin-transform-runtime']
                    }
                },
            }
        ]
    },
    plugins: [
        new CleanWebpackPlugin(),
        new HtmlWebpackPlugin({
            template: "./apps/lrr_db/lrr_db.ejs",
            filename: "lrr_db.html",
            chunks: ['lrr_db']
        }),
        new HtmlWebpackPlugin({
            template: "./apps/find_lrr/find_lrr.ejs",
            filename: "find_lrr.html",
            chunks: ['find_lrr']
        }),
        new HtmlWebpackPlugin({
            template: "./apps/blog/about.ejs",
            filename: "about.html",
            mdHTML: aboutMdHtml,
            chunks: ['']
        }),
        new HtmlWebpackPlugin({
            template: "./apps/error_page/404.ejs",
            filename: "404.html",
            chunks: []
        }),
        new HtmlWebpackPlugin({
            template: "./apps/error_page/403.ejs",
            filename: "403.html",
            chunks: []
        }),
        new HtmlWebpackPlugin({
            template: "./apps/error_page/429.ejs",
            filename: "429.html",
            chunks: []
        }),
        new CopyWebpackPlugin([
            {
                context: "./doc/images/",
                from: "*",
                to: path.resolve(DIST_PATH, 'images')
            }
        ]),
        new ZipWebpackPlugin({
            filename: "dist.zip"
        }),
    ]
};