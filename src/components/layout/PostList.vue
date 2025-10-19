<!-- PostList.vue -->
<template>
  <div>
    <!-- 文章列表 -->
    <div v-for="post in displayedPosts" :key="post.id" class="col-12 mb-3">
      <div class="card shadow-sm border-0" :style="{ backgroundColor: 'var(--app-card-bg)' }">
        <div class="card-body p-4 d-flex flex-column flex-md-row gap-4">
          <!-- 内容区域 -->
          <div class="flex-grow-1">
            <!-- 日期 -->
            <small class="meta-text mb-2 d-block" :style="{ color: 'var(--app-text-muted)' }">
              {{ formatDate(post.date) }}
            </small>

            <!-- 标题 -->
            <h5 class="post-title fw-bold mb-2" :style="{ color: 'var(--app-text)' }">
              <router-link :to="getArticlePath(post)" class="text-decoration-none" :style="{ color: 'var(--app-text)' }">
                {{ post.title }}
              </router-link>
            </h5>

            <!-- 分类 -->
            <div class="mb-2">
              <small class="meta-text" :style="{ color: 'var(--app-text-muted)' }">
                {{ formatCategory(post.category) }}
              </small>
            </div>

            <!-- 预览内容 -->
            <p class="mb-3 desc-text">
              {{ post.preview }}
            </p>

            <!-- 标签 -->
            <div class="d-flex flex-wrap gap-2">
              <span v-for="tag in post.tags" :key="tag" class="badge tag-badge fw-normal py-1 px-2 rounded-3" style="cursor: pointer" @click="goTag(tag)">
                # {{ tag }}
              </span>
            </div>
          </div>
        </div>
      </div>
    </div>

    <!-- 分页组件 -->
    <div class="row mt-4 pb-4" v-if="totalPages > 1">
      <div class="col-12">
        <nav aria-label="文章列表分页">
          <ul class="pagination justify-content-between align-items-center mb-0">
            <!-- 上一页按钮 -->
            <li class="page-item" :class="{ disabled: currentPage === 1 }">
              <button class="page-link d-flex align-items-center border-0 bg-transparent px-3 py-2" :style="{ color: 'var(--app-text-muted)' }"
                @click="prevPage" :disabled="currentPage === 1" aria-label="上一页">
                <span class="me-1">&lt;</span>
                <span>上一页</span>
              </button>
            </li>

            <!-- 页码 -->
            <li class="page-item">
              <div class="d-flex gap-2">
                <!-- 第一页 -->
                <button v-if="showFirstPage" class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === 1 }" @click="goToPage(1)">
                  1
                </button>

                <!-- 前面的省略号 -->
                <span class="d-flex align-items-center" v-if="showFirstEllipsis">...</span>

                <!-- 中间的页码 -->
                <button v-for="page in middlePages" :key="page"
                  class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === page }" @click="goToPage(page)">
                  {{ page }}
                </button>

                <!-- 后面的省略号 -->
                <span class="d-flex align-items-center" v-if="showLastEllipsis">...</span>

                <!-- 最后一页 -->
                <button v-if="showLastPage && totalPages > 1"
                  class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === totalPages }" @click="goToPage(totalPages)">
                  {{ totalPages }}
                </button>
              </div>
            </li>

            <!-- 下一页按钮 -->
            <li class="page-item" :class="{ disabled: currentPage === totalPages }">
              <button class="page-link d-flex align-items-center border-0 bg-transparent px-3 py-2" :style="{ color: 'var(--app-text-muted)' }"
                @click="nextPage" :disabled="currentPage === totalPages" aria-label="下一页">
                <span>下一页</span>
                <span class="ms-1">&gt;</span>
              </button>
            </li>
          </ul>
        </nav>
      </div>
    </div>
  </div>
</template>

<script>
import notesFlat from '@/content/notes.json'
import categoriesData from '@/content/categories.json'

export default {
  name: 'PostList',
  props: {
    docs: {
      type: Array,
      required: true,
      default: () => []
    },
    perPage: {
      type: Number,
      default: 6
    }
  },
  data() {
    return {
      currentPage: 1,
      maxVisiblePages: 5
    }
  },
  computed: {
    totalPages() {
      return Math.max(1, Math.ceil(this.docs.length / this.perPage))
    },
    displayedPosts() {
      const start = (this.currentPage - 1) * this.perPage
      const end = start + this.perPage
      return this.docs.slice(start, end)
    },
    allVisiblePages() {
      const pages = []
      const total = this.totalPages
      const current = this.currentPage
      const maxShow = this.maxVisiblePages

      if (total <= maxShow) {
        for (let i = 1; i <= total; i++) {
          pages.push(i)
        }
        return pages
      }

      pages.push(current)

      let i = 1
      while (pages.length < maxShow) {
        if (current - i >= 1 && !pages.includes(current - i)) {
          pages.push(current - i)
        }
        if (pages.length >= maxShow) break

        if (current + i <= total && !pages.includes(current + i)) {
          pages.push(current + i)
        }
        i++
      }

      if (!pages.includes(1)) {
        pages.pop()
        pages.push(1)
      }
      if (pages.length < maxShow && !pages.includes(total)) {
        pages.push(total)
      } else if (pages.length >= maxShow && !pages.includes(total)) {
        pages.pop()
        pages.push(total)
      }

      return pages.sort((a, b) => a - b)
    },
    middlePages() {
      return this.allVisiblePages.filter(page =>
        page !== 1 && page !== this.totalPages
      )
    },
    showFirstPage() {
      return this.allVisiblePages.includes(1)
    },
    showLastPage() {
      return this.allVisiblePages.includes(this.totalPages)
    },
    showFirstEllipsis() {
      return this.totalPages > this.maxVisiblePages && this.allVisiblePages[1] > 2
    },
    showLastEllipsis() {
      return this.totalPages > this.maxVisiblePages &&
        this.allVisiblePages[this.allVisiblePages.length - 2] < this.totalPages - 1
    },
    // 将 categories.json 扁平为 key -> title 映射，用于把二级分类 key 显示为中文标题
    categoryTitleMap() {
      const map = {}
      try {
        (categoriesData || []).forEach(section => {
          (section.items || []).forEach(item => {
            (item.categories || []).forEach(cat => {
              if (cat && cat.key && cat.title) {
                map[cat.key] = cat.title
              }
            })
          })
        })
      } catch (e) {
        // 静默失败，保持空映射
      }
      return map
    }
  },
  methods: {
    formatDate(dateString) {
      const options = { year: 'numeric', month: 'long', day: 'numeric' }
      return new Date(dateString).toLocaleDateString('zh-CN', options)
    },
    getArticlePath(post) {
      // 根据文章标题从notes.json中查找对应的文件路径
      let articlePath = 'notes/Programming/python/25-09-19--python-fastq.md' // 默认路径
      
      // 在扁平 notes.json 中查找匹配的文章
      const found = Array.isArray(notesFlat) ? notesFlat.find(item => item.title === post.title) : null
      if (found && found.relativePath) {
        articlePath = `notes/${found.relativePath}.md`
      }
      return `/article/${articlePath.replace(/\.md$/, '')}`
    },
    goToPage(page) {
      if (page >= 1 && page <= this.totalPages) {
        this.currentPage = page;
        const q = { ...this.$route.query, page: String(page) };
        this.$router.push({ path: this.$route.path, query: q }).catch(() => {});
        this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
      }
    },
    prevPage() {
      if (this.currentPage > 1) {
        const page = this.currentPage - 1;
        this.currentPage = page;
        const q = { ...this.$route.query, page: String(page) };
        this.$router.push({ path: this.$route.path, query: q }).catch(() => {});
        this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
      }
    },
    nextPage() {
      if (this.currentPage < this.totalPages) {
        const page = this.currentPage + 1;
        this.currentPage = page;
        const q = { ...this.$route.query, page: String(page) };
        this.$router.push({ path: this.$route.path, query: q }).catch(() => {});
        this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
      }
    },
    goTag(tag) {
      if (!tag) return;
      const q = { ...this.$route.query, tag: tag, page: '1' };
      this.$router.push({ path: '/', query: q }).catch(() => {});
      this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
    },
    clearTag() {
      const q = { ...this.$route.query };
      delete q.tag;
      q.page = '1';
      this.$router.push({ path: '/', query: q }).catch(() => {});
      this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
    },
    handleResize() {
      this.maxVisiblePages = window.innerWidth < 480 ? 3 : 5
    },
    // 显示为 "顶级 / 二级中文标题"
    formatCategory(catArr) {
      if (!Array.isArray(catArr) || catArr.length === 0) return ''
      const top = catArr[0] || ''
      const subKey = catArr[1]
      if (!subKey) return top
      const subTitle = this.categoryTitleMap[subKey] || subKey
      return `${top} / ${subTitle}`
    }
  },
  watch: {
    docs() {
      this.currentPage = 1
    },
    '$route.query.page'(newVal) {
      const p = parseInt(newVal);
      const page = Number.isFinite(p) && p >= 1 ? Math.min(p, this.totalPages) : 1;
      if (page !== this.currentPage) {
        this.currentPage = page;
        this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
      }
    }
  },
  mounted() {
    const p = parseInt(this.$route.query.page);
    this.currentPage = Number.isFinite(p) && p >= 1 ? Math.min(p, this.totalPages) : 1;
    this.handleResize();
    window.addEventListener('resize', this.handleResize);
  },
  beforeUnmount() {
    window.removeEventListener('resize', this.handleResize)
  }
}
</script>

<style scoped>
/* 标题提升 1-2 级，与文章页 h3 视觉一致 */
.post-title { font-size: 1.4rem; font-weight: 700; }

/* 描述文本与文章页正文一致 */
.desc-text { font-size: 1rem; line-height: 1.8; color: var(--app-text); }

/* 日期与分类统一为 1rem */
.meta-text { font-size: 0.95rem; }

/* 标签样式与文章页面完全一致 */
.badge { font-size: 0.85rem; font-weight: 500; }
.tag-badge { 
  color: var(--app-tag-text) !important; 
  background-color: var(--app-tag-bg) !important; 
}

.page-link {
  transition: all 0.2s ease;
  color: var(--app-pagination-link-text);
  background-color: var(--app-pagination-link-bg);
  border: none;
  outline: none;
  font-weight: 500;
  min-width: 36px;
  height: 36px;
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
}

.page-link:hover {
  background-color: var(--app-primary-bg-subtle);
  color: var(--app-primary);
}

.current-page {
  background-color: var(--app-primary) !important;
  color: var(--app-pagination-current-text) !important;
  font-weight: 600 !important;
}

.page-link:focus {
  box-shadow: none;
  outline: none;
}
</style>
