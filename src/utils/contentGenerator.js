// 自动生成内容JSON文件的工具
import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

export class ContentGenerator {
  constructor(contentDir = '../../src/content') {
    this.contentDir = path.resolve(__dirname, contentDir)
  }

  // 扫描目录结构并生成notes.json
  async generateNotesJSON() {
    const notesDir = path.join(this.contentDir, 'notes')
    const categories = []

    try {
      // 扫描一级分类（组学技术、编程语言等）
      const categoryDirs = fs.readdirSync(notesDir).filter(dir => 
        fs.statSync(path.join(notesDir, dir)).isDirectory()
      )

      for (const categoryName of categoryDirs) {
        const categoryPath = path.join(notesDir, categoryName)
        const children = []

        // 扫描二级分类
        const subcategoryDirs = fs.readdirSync(categoryPath).filter(dir => 
          fs.statSync(path.join(categoryPath, dir)).isDirectory()
        )

        for (const subcategoryName of subcategoryDirs) {
          const subcategoryPath = path.join(categoryPath, subcategoryName)
          const files = []

          // 扫描Markdown文件
          const mdFiles = fs.readdirSync(subcategoryPath).filter(file => 
            file.endsWith('.md') && !file.endsWith('-en.md')
          )

          for (const mdFile of mdFiles) {
            const filePath = path.join(subcategoryPath, mdFile)
            const content = fs.readFileSync(filePath, 'utf8')
            
            // 解析front matter获取标题
            const titleMatch = content.match(/title:\s*"([^"]+)"/)
            const title = titleMatch ? titleMatch[1] : mdFile.replace('.md', '').replace(/^\d{2}-\d{2}-\d{2}--/, '')
            
            files.push({
              title: title,
              path: `notes/${categoryName}/${subcategoryName}/${mdFile}`
            })
          }

          if (files.length > 0) {
            children.push({
              name: this.formatName(subcategoryName),
              files: files
            })
          }
        }

        if (children.length > 0) {
          categories.push({
            name: this.formatName(categoryName),
            children: children
          })
        }
      }

      return { notes: categories }
    } catch (error) {
      console.error('Error generating notes JSON:', error)
      return { notes: [] }
    }
  }

  // 生成posts.json（从notes.json提取）
  async generatePostsJSON(notesData) {
    const posts = []
    let id = 1

    for (const category of notesData.notes) {
      for (const subcategory of category.children) {
        for (const file of subcategory.files) {
          const filePath = path.join(this.contentDir, file.path)
          
          try {
            const content = fs.readFileSync(filePath, 'utf8')
            const dateMatch = content.match(/date:\s*"([^"]+)"/)
            const authorMatch = content.match(/author:\s*"([^"]+)"/)
            const tagsMatch = content.match(/tags:\s*\[([^\]]+)\]/)
            const descriptionMatch = content.match(/description:\s*"([^"]+)"/)

            const tags = tagsMatch ? 
              tagsMatch[1].split(',').map(tag => tag.trim().replace(/"/g, '')) : 
              []

            posts.push({
              id: id++,
              title: file.title,
              date: dateMatch ? dateMatch[1] : '2025-09-19',
              category: [category.name, subcategory.name],
              tags: tags,
              preview: descriptionMatch ? descriptionMatch[1] : `${file.title}的详细内容`,
              path: file.path
            })
          } catch (error) {
            console.error(`Error processing file ${filePath}:`, error)
          }
        }
      }
    }

    return posts
  }

  // 格式化名称（将文件夹名转换为友好名称）
  formatName(name) {
    const nameMap = {
      'Programming': '编程语言',
      'Bioinformatics': '生物信息学',
      'Omics': '组学技术',
      'DataScience': '数据科学',
      'python': 'Python',
      'r': 'R语言',
      'shell': 'Shell',
      'javascript': 'JavaScript',
      'alignment': '序列比对',
      'structure': '结构分析',
      'genomics': '基因组学',
      'proteomics': '蛋白质组学',
      'transcriptomics': '转录组学',
      'statistics': '统计分析',
      'machinelearning': '机器学习',
      'visualization': '数据可视化'
    }

    return nameMap[name] || name
  }

  // 生成所有JSON文件
  async generateAll() {
    const notesData = await this.generateNotesJSON()
    const postsData = await this.generatePostsJSON(notesData)

    // 写入notes.json
    fs.writeFileSync(
      path.join(this.contentDir, 'notes', 'notes.json'),
      JSON.stringify(notesData, null, 2)
    )

    // 写入posts.json
    fs.writeFileSync(
      path.join(this.contentDir, 'posts.json'),
      JSON.stringify(postsData, null, 2)
    )

    console.log('Content JSON files generated successfully!')
    return { notesData, postsData }
  }
}

// 如果直接运行此文件，则生成JSON
if (import.meta.url === `file://${process.argv[1]}`) {
  const generator = new ContentGenerator()
  generator.generateAll().catch(console.error)
}