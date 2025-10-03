<template>
  <div class="navigation-tree">
    <!-- 只显示当前文章所在的一级分类 -->
    <div v-if="currentCategoryData" class="category-section">

      <!-- 一级分类标题 -->
      <div class="category-header">
        <h3 class="category-title">{{ currentCategoryData.name }}</h3>
      </div>

      <!-- 二级分类和文章列表 -->
      <div class="category-content">
        <div v-for="subcategory in currentCategoryData.children" :key="subcategory.name" class="subcategory-section">
          <!-- 二级分类标题 -->
          <div class="subcategory-header">
            <h4 class="subcategory-title" :class="{ 'current-subcategory': isCurrentSubcategory(subcategory.name) }">
              {{ subcategory.name }}
            </h4>
          </div>
          
          <!-- 文章列表 -->
          <div class="article-list-container">
            <ul class="article-list">
              <li v-for="file in subcategory.files" :key="file.title" class="article-item">
                <router-link :to="getArticlePath(file)" class="article-link" :class="{ 'current-article': isCurrentArticle(file) }">
                  <span class="article-title">{{ file.title }}</span>
                </router-link>
              </li>
            </ul>
          </div>
        </div>
      </div>
    </div>
    
    <!-- 如果没有当前文章分类，显示提示 -->
    <div v-else class="no-category-message">
      <p class="text-muted">请选择一篇文章查看目录</p>
    </div>
  </div>
</template>

<script>
import notesData from '@/content/notes/notes.json'

export default {
  name: 'NavigationTree',
  data() {
    return {
      notesData: notesData.notes || [],
      expandedSubcategories: {},
      currentCategory: '',
      currentSubcategory: '',
      currentArticle: ''
    }
  },
  computed: {
    currentCategoryData() {
      if (!this.currentCategory) return null
      return this.notesData.find(cat => cat.name === this.currentCategory)
    }
  },
  mounted() {
    // 根据当前路由确定当前文章的分类
    this.determineCurrentArticle()
    
    // 默认所有二级分类都折叠，只有当前文章所在的二级分类展开
    if (this.currentCategoryData) {
      this.currentCategoryData.children.forEach(subcat => {
        this.expandedSubcategories[subcat.name] = true
      })
    }
  },
  watch: {
    '$route': {
      handler() {
        this.determineCurrentArticle()
        // 更新展开状态
        if (this.currentCategoryData) {
          this.currentCategoryData.children.forEach(subcat => {
            this.expandedSubcategories[subcat.name] = true
          })
        }
      },
      immediate: true,
      deep: true
    }
  },
  methods: {
    toggleSubcategory(subcategoryName) {
      // 切换指定二级分类的展开状态
      this.expandedSubcategories[subcategoryName] = !this.expandedSubcategories[subcategoryName]
    },
    getArticlePath(file) {
      // 生成完整的路由路径
      const path = file.path.replace(/\.md$/, '')
      return `/article/${path}`
    },
    isCurrentArticle(file) {
      return file.path === this.currentArticle
    },
    isCurrentSubcategory(subcategoryName) {
      return subcategoryName === this.currentSubcategory
    },
    determineCurrentArticle() {
      // 根据当前路由参数确定当前文章
      const route = this.$route
      // The path param is an array of segments, join them
      const path = route.params.path ? route.params.path.join('/') : ''
      
      if (path) {
        // 查找匹配的文章路径
        for (const category of this.notesData) {
          for (const subcategory of category.children) {
            const foundFile = subcategory.files.find(file => 
              file.path.replace(/\.md$/, '') === path
            )
            if (foundFile) {
              this.currentCategory = category.name
              this.currentSubcategory = subcategory.name
              this.currentArticle = foundFile.path
              return
            }
          }
        }
      }
      
      // 如果没有路由参数或未找到匹配，默认显示FASTQ文章
      const programmingCategory = this.notesData.find(cat => cat.name === "编程语言")
      if (programmingCategory) {
        const pythonSubcat = programmingCategory.children.find(sub => sub.name === "Python")
        if (pythonSubcat) {
          const fastqFile = pythonSubcat.files.find(file => 
            file.path.includes('python-fastq')
          )
          if (fastqFile) {
            this.currentCategory = "编程语言"
            this.currentSubcategory = "Python"
            this.currentArticle = fastqFile.path
            return
          }
        }
      }
      
      // 如果连FASTQ都找不到，使用第一个可用的文章
      if (this.notesData.length > 0 && this.notesData[0].children.length > 0) {
        const firstSubcat = this.notesData[0].children[0]
        if (firstSubcat.files.length > 0) {
          this.currentCategory = this.notesData[0].name
          this.currentSubcategory = firstSubcat.name
          this.currentArticle = firstSubcat.files[0].path
        }
      }
    }
  }
}
</script>

<style scoped>
.navigation-tree {
  font-family: inherit;
}

.category-header {
  margin-bottom: 0.5rem;
}

.category-title {
  font-size: 1rem;
  font-weight: 700;
  color: var(--bs-secondary-color);
  margin: 0;
}

.subcategory-section {
  margin-bottom: 0rem;
}

.subcategory-header {
  cursor: default;
  padding: 0.25rem 0;
}

.subcategory-title {
  font-size: 0.95rem;
  font-weight: 700;
  color: var(--bs-body-color);
  margin: 0;
  display: flex;
  justify-content: space-between;
  align-items: center;
  transition: color 0.2s ease;
}

.article-list-container {
  padding: 0.1rem 0 0 0;
  margin-left: 0.5rem;
}

.article-list {
  list-style: none;
  padding: 0;
  margin: 0;
}

.article-item {
  margin-bottom: 0.25rem;
}

.article-link {
  display: block;
  padding: 0.4rem 0.75rem;
  color: var(--bs-body-color);
  text-decoration: none;
  border-radius: 0.25rem;
  transition: all 0.2s ease;
  font-size: 0.9rem;
  margin-left: 0.25rem;
}



.article-link.current-article {
  color: var(--bs-primary);
  font-weight: 700;
}

.article-title {
  margin: 0;
  line-height: 1.4;
}

.no-category-message {
  text-align: center;
  padding: 2rem 1rem;
  color: var(--bs-secondary-color);
  font-size: 0.9rem;
}
</style>
