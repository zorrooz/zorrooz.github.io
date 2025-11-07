<!-- Category.vue -->
<template>
  <div class="container py-4 px-3 view-container category-view">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8 typography-body">
        <!-- 页面标题 -->
        <div class="text-center mb-4">
          <h1 class="article-title mb-3">{{ pageTitle }}</h1>
        </div>

        <!-- 分类模块 -->
        <div class="d-flex flex-column" style="gap: 2rem;">
          <section v-for="(category, index) in categoryList" :key="index">
            <!-- 一级分类标题：浅灰下划线 -->
            <h2 class="h4 fw-semibold pb-2 mb-3 heading-underline" style="color: var(--app-text-emphasis);">
              {{ category.title }}
            </h2>

            <!-- 卡片网格 -->
            <div class="row g-3">
              <div v-for="(item, idx) in category.items" :key="idx" class="col-12 col-md-6 col-lg-4">
                <div class="card h-100 shadow-sm border-0" :style="{ backgroundColor: 'var(--app-card-bg)' }">
                  <div class="card-body p-4 d-flex flex-column">
                    <!-- 标题 -->
                    <h3 class="h5 fw-semibold mb-1" :style="{ color: 'var(--app-text-emphasis)' }">
                      {{ item.title || item.name }}
                    </h3>
                    <!-- 信息行：日期紧跟图标（项目/课题），笔记仅日期 -->
                    <div v-if="getLatestDate(item)" class="d-inline-flex align-items-center gap-2 mb-2">
                      <i class="fas fa-calendar-alt meta-icon"></i>
                      <small class="meta-text" :style="{ color: 'var(--app-text-muted)' }">{{ formatMonth(getLatestDate(item)) }}</small>
                      <!-- 项目：GitHub 图标紧跟日期 -->
                      <a
                        v-if="category.title === projectsTitle && item.github"
                        :href="normalizeUrl(item.github)"
                        target="_blank"
                        rel="noopener noreferrer"
                        class="icon-link"
                        title="GitHub"
                        @click.stop
                      >
                        <i class="fab fa-github"></i>
                      </a>
                      <!-- 课题：链接图标紧跟日期 -->
                      <a
                        v-else-if="category.title === topicsTitle && item.doi"
                        :href="normalizeUrl(item.doi)"
                        target="_blank"
                        rel="noopener noreferrer"
                        class="icon-link"
                        title="链接"
                        @click.stop
                      >
                        <i class="fas fa-link"></i>
                      </a>
                    </div>

                    <!-- 动作图标并入标题下方的信息行，此处移除 -->

                    <!-- 简介 -->
                    <p class="mb-2 flex-grow-1 desc-text" :style="{ color: 'var(--app-text-secondary)' }">
                      {{ item.desc }}
                    </p>

                    <!-- 标签（项目与课题）：与 PostList 一致的轻胶囊 -->
                    <div v-if="category.title !== notesTitle && Array.isArray(item.tags) && item.tags.length" class="d-flex flex-wrap gap-2 mb-2">
                      <span v-for="(tag, tIdx) in item.tags" :key="tIdx" class="badge tag-badge fw-normal py-1 px-2 rounded-3">
                        # {{ tag }}
                      </span>
                    </div>

                    <!-- 指标区（笔记）：完全参考 ProfileCard 统计条风格 -->
                    <div v-if="category.title === t('notes')" class="row g-0 text-center stats-row mb-2 w-100" :style="{ 'border-color': 'var(--app-border)' }">
                      <div v-if="Array.isArray(item.tags) && item.tags.length" class="col border-end">
                        <div class="fw-bold stat-num">{{ item.tags.length }}</div>
                        <div class="stat-label" :style="{ color: 'var(--app-text-muted)' }">{{ t('tags') }}</div>
                      </div>
                      <div v-if="item.stats && item.stats.postsCount" class="col border-end">
                        <div class="fw-bold stat-num">{{ item.stats.postsCount }}</div>
                        <div class="stat-label" :style="{ color: 'var(--app-text-muted)' }">{{ articlesText }}</div>
                      </div>
                      <div v-if="item.stats && item.stats.totalWords" class="col">
                        <div class="fw-bold stat-num">{{ formatWords(item.stats.totalWords) }}</div>
                        <div class="stat-label" :style="{ color: 'var(--app-text-muted)' }">{{ wordsText }}</div>
                      </div>
                    </div>

                    <!-- 元信息区：克制风格，图标与少量文字 -->
                    <div class="d-flex flex-wrap align-items-center gap-2 mb-2 meta-row">

                    </div>

                    <!-- 分隔线：更浅的灰色 -->
                    <hr class="my-2" />

                    <!-- 查看更多 -->
                    <div class="text-end">
                      <span class="fw-medium d-inline-flex align-items-center gap-1 cursor-pointer see-more-text"
                        :style="{ color: 'var(--app-primary)', transition: 'transform 0.2s ease' }" @click="() => handleSeeMore(item)"
                        @mouseenter="e => e.target.style.transform = 'translateX(2px)'"
                        @mouseleave="e => e.target.style.transform = 'translateX(0)'">
                        {{ seeMoreText }}
                        <i class="bi bi-arrow-right"></i>
                      </span>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          </section>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import { useI18n } from 'vue-i18n';
import { loadCategories } from '@/utils/contentLoader';


export default {
  name: 'CategoryView',
  setup() {
    const { t, locale } = useI18n();
    return { t, locale };
  },
  data() {
    return {
      categoryList: []
    };
  },
  computed: {
    pageTitle() {
      return this.t('categories');
    },
    notesTitle() {
      return this.t('notes');
    },
    projectsTitle() {
      return this.t('projects');
    },
    topicsTitle() {
      return this.t('topics');
    },
    tagsText() {
      return this.t('tags');
    },
    articlesText() {
      return this.t('articles');
    },
    wordsText() {
      return this.t('words');
    },
    seeMoreText() {
      return this.t('seeMore');
    }
  },
  async created() {
    await this.loadCategoryData();
  },
  
  watch: {
    locale() {
      this.loadCategoryData();
    }
  },
  
  methods: {
    async loadCategoryData() {
      try {
        const categoryData = await loadCategories();
        this.categoryList = Array.isArray(categoryData) ? categoryData : (categoryData?.categoryList || []);
      } catch (error) {
        console.error('Failed to load category data:', error);
        this.categoryList = [];
      }
    },
    getLatestDate(item) {
      return (item && item.stats && item.stats.latestDate) || '';
    },
    normalizeUrl(u) {
      if (typeof u !== 'string' || !u.trim()) return '';
      if (/^https?:\/\//i.test(u)) return u;
      return 'https://' + u.replace(/^\/+/, '');
    },
    formatMonth(dateStr) {
      // 仅显示到月份：YYYY-MM
      const d = new Date(dateStr);
      if (isNaN(d.getTime())) return dateStr;
      const y = d.getFullYear();
      const m = String(d.getMonth() + 1).padStart(2, '0');
      return `${y}-${m}`;
    },
    formatWords(n) {
      // 简单缩写：>=10k 用 k 显示
      const num = Number(n);
      if (!Number.isFinite(num)) return n;
      if (num >= 10000) return Math.round(num / 1000) + 'k';
      return String(num);
    },
    handleSeeMore(item) {
      const root = item?.root;
      if (root) {
        this.$router.push(root).catch(err => {
          if (err.name !== 'NavigationDuplicated' && !err.toString().includes('Navigation cancelled')) {
            console.error('Navigation error:', err);
          }
        });
        return;
      }

      // 无 root 时，若有 doi 外链则打开（自动补协议）
      const doi = this.normalizeUrl(item?.doi);
      if (doi) {
        window.open(doi, '_blank', 'noopener,noreferrer');
      }
    }
  }
};
</script>

<style scoped>
.typography-body { font-size: 1.125rem; line-height: 1.8; }
.typography-body p { margin-bottom: 0.75rem; }
.heading-underline { border-bottom: 1px solid var(--app-border); }

.desc-text { font-size: 1rem; line-height: 1.6; }
.see-more-text { font-size: 0.9rem; }

/* 动作图标样式：与 Header 一致的扁平风格 */
.icon-link { 
  color: var(--app-icon-color); 
  font-size: 1.1rem; 
  display: inline-flex; 
  align-items: center; 
  justify-content: center;
  width: 24px;
  height: 24px;
  text-decoration: none; 
}
.icon-link:hover { color: var(--app-icon-hover-color); text-decoration: none; }
/* 统一图标尺寸与点击区域，确保 GitHub 与链接一致 */
.icon-link i { font-size: 1rem; line-height: 1; }

/* 元信息区样式 */
.meta-row { row-gap: 0.25rem; }
.meta-text { font-size: 0.95rem; }
.meta-icon { font-size: 0.95rem; color: var(--app-text-muted); }

/* 指标区（笔记）：ProfileCard 统计条风格（更紧凑） */
.stats-row { }
.stats-row .col { padding: 4px 0; }
.stats-row .col:not(:last-child) { border-right: 1px solid var(--app-border); }
.stat-num { font-size: 1rem; line-height: 1.2; font-weight: 700; color: var(--app-stat-num-color); margin-bottom: 0; }
.stat-label { font-size: 0.85rem; line-height: 1.2; color: var(--app-text-muted); margin-top: 0; }

/* 标签与 TagCloud 保持一致 */
.badge { font-size: 0.85rem; font-weight: 500; }
.tag-badge { 
  color: var(--app-tag-text) !important; 
  background-color: var(--app-tag-bg) !important; 
  font-size: 0.85rem !important;
}

.card { border-radius: 0.75rem; }
.card-body { padding: 1.5rem !important; }

</style>
