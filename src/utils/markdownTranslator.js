import OpenAI from "openai";
import fs from "fs/promises";
import path from "path";

const openai = new OpenAI({
  baseURL: 'https://api.deepseek.com',
  apiKey: '' // 🔑 替换为你的实际 API Key
});

async function translate(text) {
  try {
    const completion = await openai.chat.completions.create({
      messages: [
        {
          role: "system",
          content: "你是一个专业翻译器。请将以下内容翻译为英文，要求：\n1. 严格保持原始格式（包括 Markdown 语法、换行、缩进、标点等）\n2. 不添加任何解释、注释或额外内容\n3. 只输出翻译结果"
        },
        {
          role: "user",
          content: text
        }
      ],
      model: "deepseek-chat",
      temperature: 0.3
    });

    return completion.choices[0].message.content.trim();
  } catch (error) {
    console.error("翻译时发生错误:", error.message);
    throw error;
  }
}

async function translateFile(inputFilePath) {
  try {
    // 读取文件内容
    const content = await fs.readFile(inputFilePath, "utf-8");

    // 翻译
    console.log(`正在翻译文件: ${inputFilePath}`);
    const translated = await translate(content);

    // 生成输出路径: example.md → example-en.md
    const ext = path.extname(inputFilePath);
    const basename = path.basename(inputFilePath, ext);
    const outputPath = path.join(
      path.dirname(inputFilePath),
      `${basename}-en${ext}`
    );

    // 写入翻译后文件
    await fs.writeFile(outputPath, translated, "utf-8");
    console.log(`✅ 翻译完成，已保
    存至: ${outputPath}`);
  } catch (error) {
    console.error("处理文件时出错:", error.message);
  }
}

// 🧑‍💻 在这里手动指定你要翻译的文件路径
async function main() {
  const filePath = "../content/article.md"; // ← ← ← 你在这里改路径！

  if (!filePath.endsWith(".md")) {
    console.log("⚠️  请确保路径指向 .md 文件");
    return;
  }

  await translateFile(filePath);
}

main();
