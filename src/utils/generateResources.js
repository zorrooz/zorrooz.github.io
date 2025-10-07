import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import yaml from 'js-yaml';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentSrcDir = path.join(__dirname, '../content-src');
const contentOutputDir = path.join(__dirname, '../content');

function ensureDirectoryExistence(filePath) {
  const dirname = path.dirname(filePath);
  if (!fs.existsSync(dirname)) {
    fs.mkdirSync(dirname, { recursive: true });
  }
}

// 数据结构（严格解析当前结构，仅做安全归一化）：
// Array<{
//   title: string,
//   children: Array<{
//     title: string,
//     items: Array<{
//       name: string,
//       url: string,
//       desc?: string
//     }>
//   }>
// }>
function normalize(raw) {
  const list = Array.isArray(raw) ? raw : [];
  return list.map(cat => {
    const title = typeof cat?.title === 'string' ? cat.title : '';
    const children = Array.isArray(cat?.children) ? cat.children : [];
    const normChildren = children.map(sub => {
      const st = typeof sub?.title === 'string' ? sub.title : '';
      const items = Array.isArray(sub?.items) ? sub.items : [];
      const normItems = items.map(it => ({
        name: typeof it?.name === 'string' ? it.name : '',
        url: typeof it?.url === 'string' ? it.url : '',
        desc: typeof it?.desc === 'string' ? it.desc : ''
      }));
      return { title: st, items: normItems };
    });
    return { title, children: normChildren };
  });
}

function generateResourcesJson() {
  const yamlPath = path.join(contentSrcDir, 'resources.yaml');
  const outputPath = path.join(contentOutputDir, 'resources.json');

  let raw = [];
  if (fs.existsSync(yamlPath)) {
    try {
      const yamlContent = fs.readFileSync(yamlPath, 'utf-8');
      const parsed = yaml.load(yamlContent);
      raw = parsed || [];
    } catch (e) {
      console.warn(`Warn: failed to parse YAML at ${yamlPath}, using empty array.`, e?.message || e);
      raw = [];
    }
  } else {
    console.warn(`Warn: source YAML not found at ${yamlPath}, using empty array.`);
  }

  const result = normalize(raw);

  try {
    ensureDirectoryExistence(outputPath);
    fs.writeFileSync(outputPath, JSON.stringify(result, null, 2), 'utf-8');
    console.log(`Successfully generated: ${outputPath}`);
  } catch (error) {
    console.error(`Failed to generate ${outputPath}:`, error);
  }
}

function main() {
  console.log('Starting resources.json generation script...');
  generateResourcesJson();
  console.log('resources.json generation complete.');
}

main();

export { generateResourcesJson };